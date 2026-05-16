// ncurses front-end for GONE2's ProgressStatus.
//
// Owns a background thread that polls Snapshot() at ~10 Hz and
// redraws the screen.  The chart uses braille line drawing (8× the
// resolution of normal cells), progress bars use Unicode 1/8th
// blocks for smooth partial fills, and colour pairs distinguish
// phases, the current task, axes, and the chart line.
//
// Linked only into the gone-ncurses target via -DGONE_NCURSES_TUI.
// Other build targets compile this file as an empty translation unit
// so the wildcard source list stays simple.

#include "ncurses_tui.hpp"

#ifdef GONE_NCURSES_TUI

#include <atomic>
#include <algorithm>
#include <chrono>
#include <clocale>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <ncurses.h>
#include <pthread.h>
#include <signal.h>
#include <unistd.h>

#include "progress.hpp"

namespace {

std::thread g_thread;
std::atomic<bool> g_running{false};
std::atomic<bool> g_stop{false};
// Set by the main thread once the GA has finished so the renderer
// switches its bottom hint from "abort" to "exit".
std::atomic<bool> g_finished{false};
// Set by the SIGINT handler installed by WaitForSigint.
volatile std::sig_atomic_t g_sigintReceived = 0;
ProgressStatus* g_progress = nullptr;
std::mutex g_screenMu;
bool g_hasColors = false;
// std::cerr capture: while ncurses owns the screen we route stderr
// into a stringstream so error messages don't get scribbled over by
// ncurses refreshes. StopNcursesTui restores the buffer and replays
// the captured text after endwin().
std::stringstream* g_cerrCapture = nullptr;
std::streambuf*    g_cerrOrigBuf = nullptr;

constexpr int kLeftMargin = 2;

// ---- Colour pair identifiers ---------------------------------------------
enum ColorId {
  CP_TITLE      = 1,
  CP_BORDER     = 2,
  CP_DONE       = 3,
  CP_CURRENT    = 4,
  CP_PENDING    = 5,
  CP_BAR_FILL   = 6,
  CP_BAR_EMPTY  = 7,
  CP_AXIS       = 8,
  CP_CURVE      = 9,
  CP_REFERENCE  = 13,  // ground-truth Ne overlay (-f flag)
  CP_STATUS     = 10,
  CP_HINT       = 11,
  CP_AXIS_LABEL = 12,
};

void InitColors() {
  if (!has_colors()) return;
  start_color();
  use_default_colors();
  init_pair(CP_TITLE,      COLOR_CYAN,    -1);
  init_pair(CP_BORDER,     COLOR_BLUE,    -1);
  init_pair(CP_DONE,       COLOR_GREEN,   -1);
  init_pair(CP_CURRENT,    COLOR_YELLOW,  -1);
  init_pair(CP_PENDING,    COLOR_WHITE,   -1);
  init_pair(CP_BAR_FILL,   COLOR_GREEN,   -1);
  init_pair(CP_BAR_EMPTY,  COLOR_WHITE,   -1);
  init_pair(CP_AXIS,       COLOR_BLUE,    -1);
  init_pair(CP_CURVE,      COLOR_CYAN,    -1);
  init_pair(CP_STATUS,     COLOR_WHITE,   -1);
  init_pair(CP_HINT,       COLOR_BLACK,   -1);
  init_pair(CP_AXIS_LABEL, COLOR_WHITE,   -1);
  init_pair(CP_REFERENCE,  COLOR_RED,     -1);
  g_hasColors = true;
}

struct Attr {
  int attr;
  Attr(int a) : attr(a) { attron(attr); }
  ~Attr() { attroff(attr); }
};

int Color(ColorId id) { return g_hasColors ? COLOR_PAIR(id) : 0; }

// ---- Helpers --------------------------------------------------------------

void DrawBox(int top, int left, int height, int width) {
  if (width < 2 || height < 2) return;
  Attr a(Color(CP_BORDER));
  mvaddstr(top, left, "╭");
  for (int x = 1; x < width - 1; ++x) mvaddstr(top, left + x, "─");
  mvaddstr(top, left + width - 1, "╮");
  for (int y = 1; y < height - 1; ++y) {
    mvaddstr(top + y, left, "│");
    mvaddstr(top + y, left + width - 1, "│");
  }
  mvaddstr(top + height - 1, left, "╰");
  for (int x = 1; x < width - 1; ++x) mvaddstr(top + height - 1, left + x, "─");
  mvaddstr(top + height - 1, left + width - 1, "╯");
}

void DrawCenteredTitle(int row, int width, const char* text) {
  const int len = static_cast<int>(std::strlen(text));
  const int trimmed = std::min(len, width);
  const int col = std::max(0, (width - trimmed) / 2);
  Attr a(Color(CP_TITLE) | A_BOLD);
  mvaddnstr(row, col, text, trimmed);
}

// Eighths blocks ▏▎▍▌▋▊▉ — 1/8 .. 7/8 width fills, with █ as 8/8.
static const char* kEighthBlocks[9] = {
  " ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█",
};

// Draw a progress bar in [col, col+width) using 1/8-block resolution.
void DrawProgressBar(int row, int col, int width, double progress) {
  if (width < 1) return;
  const double clamped = std::max(0.0, std::min(1.0, progress));
  const int totalEighths = width * 8;
  const int filledEighths =
      std::max(0, std::min(totalEighths,
                           static_cast<int>(std::round(clamped * totalEighths))));
  const int fullCells = filledEighths / 8;
  const int rem = filledEighths % 8;

  // Full cells.
  {
    Attr a(Color(CP_BAR_FILL));
    for (int i = 0; i < fullCells; ++i) mvaddstr(row, col + i, "█");
    if (rem > 0 && fullCells < width) {
      mvaddstr(row, col + fullCells, kEighthBlocks[rem]);
    }
  }
  // Empty cells.
  {
    Attr a(Color(CP_BAR_EMPTY) | A_DIM);
    const int emptyStart = fullCells + (rem > 0 ? 1 : 0);
    for (int i = emptyStart; i < width; ++i) mvaddstr(row, col + i, "·");
  }
}

void DrawPctRight(int row, int col, double pct) {
  char buf[16];
  std::snprintf(buf, sizeof(buf), "%5.1f%%", pct * 100.0);
  Attr a(Color(CP_PENDING));
  mvaddstr(row, col, buf);
}

// ---- Phase list -----------------------------------------------------------

void DrawPhases(int top, int left, int width,
                const std::vector<ProgressPhaseView>& phases) {
  int row = top;
  const int barCol  = left + 32;
  const int pctCol  = std::max(barCol + 12, width - 9);
  const int barWidth = std::max(8, pctCol - barCol - 1);

  for (const ProgressPhaseView& phase : phases) {
    ColorId cid = CP_PENDING;
    const char* glyph = "○";
    int extraAttr = A_DIM;
    if (phase.state == ProgressPhaseState::Done) {
      cid = CP_DONE; glyph = "✔"; extraAttr = 0;
    } else if (phase.state == ProgressPhaseState::Current) {
      cid = CP_CURRENT; glyph = "▶"; extraAttr = A_BOLD;
    }
    {
      Attr a(Color(cid) | extraAttr);
      mvaddstr(row, left, glyph);
      mvaddch(row, left + 1, ' ');
      const int nameMax = barCol - left - 3;
      mvaddnstr(row, left + 2, phase.name.c_str(),
                std::min<int>(phase.name.size(), nameMax));
    }
    if (phase.state != ProgressPhaseState::Pending) {
      DrawProgressBar(row, barCol, barWidth, phase.progress);
      DrawPctRight(row, pctCol, phase.progress);
    }
    ++row;
  }
}

// ---- Braille line chart ---------------------------------------------------

// Bit positions for the standard Unicode braille pattern, indexed by
// [row 0..3][col 0..1].
constexpr unsigned char kBrailleBit[4][2] = {
  {0x01, 0x08},
  {0x02, 0x10},
  {0x04, 0x20},
  {0x40, 0x80},
};

// Encode a braille bitmask into a 4-byte (UTF-8 + NUL) buffer.
void EncodeBraille(unsigned char mask, char buf[4]) {
  const unsigned int cp = 0x2800u + mask;
  buf[0] = 0xE2;
  buf[1] = 0xA0 | ((cp >> 6) & 0x3F);
  buf[2] = 0x80 | (cp & 0x3F);
  buf[3] = 0;
}

void DrawNeChart(int top, int left, int height, int width,
                 const NeEstimateView& ne, const BestScoreView& bestScore) {
  if (width < 32 || height < 8) return;

  DrawBox(top, left, height, width);
  {
    Attr a(Color(CP_TITLE) | A_BOLD);
    mvaddstr(top, left + 2, "─ Ne estimate (live) ");
  }
  if (!ne.unitLabel.empty()) {
    const int col = std::max<int>(left + 4,
                                  left + width - 4 - (int)ne.unitLabel.size());
    Attr a(Color(CP_AXIS_LABEL) | A_DIM);
    mvaddstr(top, col - 2, "─ ");
    mvaddstr(top, col, ne.unitLabel.c_str());
    mvaddstr(top, col + ne.unitLabel.size(), " ─");
  }

  // Inside-box plot area, accounting for box borders and Y-axis labels.
  // Layout:
  //   col: [border][space][y-tick label 6 cols][axis │][plot ...][space][border]
  //   row: [title row][gap][plot rows ...][axis ─][x-tick labels][gap][border]
  const int yLabelWidth = 7;
  const int plotLeft = left + 2 + yLabelWidth + 1;       // first plot col
  const int plotRight = left + width - 3;                // last plot col
  const int plotCols = plotRight - plotLeft + 1;
  const int plotTop = top + 2;
  const int plotBottom = top + height - 4;
  const int plotRows = plotBottom - plotTop + 1;
  if (plotCols < 10 || plotRows < 4) return;

  // Pinned score box: a small rectangle in the chart's top-right corner
  // that shows the running best GA score. Sized to fit at most ~30 plot
  // cols (~half a wide terminal) so the plot keeps its full vertical
  // axis labels and most of its horizontal extent. Skipped on narrow
  // charts (would crowd the plot) or short charts (no vertical room).
  constexpr int kBoxWidth  = 28;
  constexpr int kBoxHeight = 4;
  const bool showScoreBox = bestScore.hasValue &&
                            plotCols >= kBoxWidth + 24 &&
                            plotRows >= kBoxHeight + 4;
  const int boxRight  = left + width - 2;
  const int boxLeft   = boxRight - kBoxWidth + 1;
  const int boxTop    = top + 1;
  const int boxBottom = boxTop + kBoxHeight - 1;

  if (!ne.hasValues || ne.values.empty()) {
    Attr a(Color(CP_STATUS) | A_DIM);
    mvaddstr(top + height / 2, left + (width - 36) / 2,
             "waiting for the first GA round…");
    return;
  }

  // Y range on log10 scale. Only the first kMaxGen=150 generations
  // are charted, matching the published Ne output file.
  constexpr int kMaxGen = 150;
  std::vector<double> v;
  v.reserve(std::min<std::size_t>(ne.values.size(), kMaxGen));
  for (std::size_t i = 0;
       i < ne.values.size() && i < static_cast<std::size_t>(kMaxGen); ++i) {
    v.push_back(std::max(1.0, ne.values[i]));
  }
  // Reference (ground-truth) Ne overlay from -f flag, same gen cap.
  std::vector<double> vRef;
  vRef.reserve(std::min<std::size_t>(ne.referenceValues.size(), kMaxGen));
  for (std::size_t i = 0;
       i < ne.referenceValues.size() &&
       i < static_cast<std::size_t>(kMaxGen); ++i) {
    vRef.push_back(std::max(1.0, ne.referenceValues[i]));
  }
  if (v.empty()) return;
  double mn = *std::min_element(v.begin(), v.end());
  double mx = *std::max_element(v.begin(), v.end());
  if (!vRef.empty()) {
    mn = std::min(mn, *std::min_element(vRef.begin(), vRef.end()));
    mx = std::max(mx, *std::max_element(vRef.begin(), vRef.end()));
  }
  // Snap the Y range to integer log10 decade boundaries so the axis
  // is stable across frames and tick labels land on round powers of
  // ten (… 1, 10, 100, 1000, 10000 …). Aim for a 2-decade window
  // with the top ≈ 10× the curve's max — i.e. the data's max sits
  // near the middle of the visible range, rather than at whichever
  // edge the previous floor/ceil rule happened to round to. The
  // +0.5 fudge centres max ≈ 10ⁿ exactly at 50 % and pushes
  // off-decade values (e.g. 200 or 999) toward the middle instead
  // of pinning them to the bottom of [100, 10000].
  int decadeMax = static_cast<int>(std::ceil(std::log10(mx) + 0.5));
  int decadeMin = decadeMax - 2;
  // Extend the bottom if the curve's min would clip out of view.
  const int dataMinDecade = static_cast<int>(std::floor(std::log10(mn)));
  if (dataMinDecade < decadeMin) decadeMin = dataMinDecade;
  // Clamp to the displayable range; if the top hits the ceiling,
  // pull the bottom up so we still show at least 2 decades.
  if (decadeMin < 0) decadeMin = 0;
  if (decadeMax > 7) decadeMax = 7;
  if (decadeMax - decadeMin < 2) {
    decadeMin = std::max(0, decadeMax - 2);
  }
  const double lmn = static_cast<double>(decadeMin);
  const double lmx = static_cast<double>(decadeMax);

  // Vertical axis with 4 tick labels.
  {
    Attr a(Color(CP_AXIS));
    for (int r = 0; r < plotRows; ++r) {
      mvaddstr(plotTop + r, plotLeft - 1, "│");
    }
  }
  // Tick labels at every integer log10 decade between decadeMin and
  // decadeMax (e.g. 1, 10, 100, 1k, 10k, …). Skipping minor ticks
  // keeps the labels stable when the data jitters by <1 decade.
  {
    Attr a(Color(CP_AXIS_LABEL) | A_DIM);
    for (int d = decadeMin; d <= decadeMax; ++d) {
      const double frac = (lmx == lmn) ? 1.0 : (d - lmn) / (lmx - lmn);
      const int row =
          plotTop + static_cast<int>(std::round((1.0 - frac) * (plotRows - 1)));
      if (row < plotTop || row > plotBottom) continue;
      const double linv = std::pow(10.0, static_cast<double>(d));
      char buf[16];
      if (d >= 4)        std::snprintf(buf, sizeof(buf), "%5.0e", linv);
      else               std::snprintf(buf, sizeof(buf), "%5.0f", linv);
      mvaddstr(row, plotLeft - 1 - yLabelWidth + 2, buf);
      mvaddstr(row, plotLeft - 1, "┤");
    }
  }
  // Gen-1 (most recent) Ne marker on the left axis. Painted on top
  // of any decade label that lands at the same row so the user can
  // always read the curve's starting value, in the curve's own
  // colour for an at-a-glance correspondence with the chart line.
  {
    const double gen1 = v[0];
    const double frac =
        (lmx == lmn) ? 1.0 : (std::log10(gen1) - lmn) / (lmx - lmn);
    const int row =
        plotTop + static_cast<int>(std::round((1.0 - frac) * (plotRows - 1)));
    if (row >= plotTop && row <= plotBottom) {
      char buf[16];
      if (gen1 >= 1e4) std::snprintf(buf, sizeof(buf), "%5.0e", gen1);
      else             std::snprintf(buf, sizeof(buf), "%5.0f", gen1);
      {
        Attr a(Color(CP_CURVE) | A_BOLD);
        mvaddstr(row, plotLeft - 1 - yLabelWidth + 2, buf);
        mvaddstr(row, plotLeft - 1, "►");
      }
    }
  }

  // Allocate pixel-cell bitmap. Each text cell holds a 2x4 braille
  // sub-grid, so resolution is (plotCols*2) × (plotRows*4) pixels.
  std::vector<unsigned char> cells(plotCols * plotRows, 0);
  const int pixW = plotCols * 2;
  const int pixH = plotRows * 4;
  const int N = static_cast<int>(v.size());

  // Sample the log10(Ne) curve at any sub-pixel x position by linearly
  // interpolating between adjacent data points.
  auto sampleYAt = [&](double xPix) -> int {
    if (N <= 1) {
      const double frac = (std::log10(v[0]) - lmn) /
                          std::max(1e-9, (lmx - lmn));
      return std::max(0, std::min(pixH - 1,
                                  (int)std::round((1.0 - frac) * (pixH - 1))));
    }
    double t = xPix * (N - 1.0) / std::max(1, pixW - 1);
    if (t < 0) t = 0;
    if (t > N - 1) t = N - 1;
    const int i = std::min(N - 2, std::max(0, (int)t));
    const double frac = t - i;
    const double lv = std::log10(v[i]) +
                      frac * (std::log10(v[i + 1]) - std::log10(v[i]));
    const double yFrac = (lv - lmn) / std::max(1e-9, (lmx - lmn));
    int y = (int)std::round((1.0 - yFrac) * (pixH - 1));
    if (y < 0) y = 0;
    if (y >= pixH) y = pixH - 1;
    return y;
  };

  // Translate the pinned score box into plot-cell coordinates so the
  // braille writers can skip any cell that would otherwise paint
  // pixels behind it. Cells fully inside [boxRow0..boxRow1] ×
  // [boxCol0..boxCol1] are off-limits for the curve.
  const int boxCol0 = showScoreBox ? boxLeft  - plotLeft : 0;
  const int boxCol1 = showScoreBox ? boxRight - plotLeft : -1;
  const int boxRow0 = showScoreBox ? boxTop    - plotTop : 0;
  const int boxRow1 = showScoreBox ? boxBottom - plotTop : -1;
  auto inScoreBox = [&](int textRow, int textCol) {
    return showScoreBox &&
           textRow >= boxRow0 && textRow <= boxRow1 &&
           textCol >= boxCol0 && textCol <= boxCol1;
  };

  auto setPix = [&](int x, int y) {
    if (x < 0 || x >= pixW || y < 0 || y >= pixH) return;
    const int textRow = y / 4;
    const int textCol = x / 2;
    if (inScoreBox(textRow, textCol)) return;
    cells[textRow * plotCols + textCol] |= kBrailleBit[y % 4][x % 2];
  };

  // For each pixel column, find the y-range the line occupies between
  // the previous and next sub-pixel x positions, and fill every pixel
  // in that range. The result is a solid, always-connected line that
  // never breaks across steep slopes.
  for (int xp = 0; xp < pixW; ++xp) {
    const double xLeft  = std::max(0.0, xp - 0.5);
    const double xRight = std::min((double)(pixW - 1), xp + 0.5);
    const int yL = sampleYAt(xLeft);
    const int yC = sampleYAt(xp);
    const int yR = sampleYAt(xRight);
    const int yMin = std::min({yL, yC, yR});
    const int yMax = std::max({yL, yC, yR});
    for (int y = yMin; y <= yMax; ++y) setPix(xp, y);
  }

  // Reference (ground-truth) Ne overlay — same sampling math against
  // its own pixel buffer, drawn FIRST in red so the live cyan curve
  // overdraws on top where they coincide. Red shines through wherever
  // the GA hasn't caught up to truth.
  if (!vRef.empty()) {
    std::vector<unsigned char> refCells(plotCols * plotRows, 0);
    const int N_ref = static_cast<int>(vRef.size());
    auto sampleRefYAt = [&](double xPix) -> int {
      if (N_ref <= 1) {
        const double frac = (std::log10(vRef[0]) - lmn) /
                            std::max(1e-9, (lmx - lmn));
        return std::max(0, std::min(pixH - 1,
                                    (int)std::round((1.0 - frac) *
                                                    (pixH - 1))));
      }
      double t = xPix * (N_ref - 1.0) / std::max(1, pixW - 1);
      if (t < 0) t = 0;
      if (t > N_ref - 1) t = N_ref - 1;
      const int i = std::min(N_ref - 2, std::max(0, (int)t));
      const double frac = t - i;
      const double lv = std::log10(vRef[i]) +
                        frac * (std::log10(vRef[i + 1]) - std::log10(vRef[i]));
      const double yFrac = (lv - lmn) / std::max(1e-9, (lmx - lmn));
      int y = (int)std::round((1.0 - yFrac) * (pixH - 1));
      if (y < 0) y = 0;
      if (y >= pixH) y = pixH - 1;
      return y;
    };
    auto setRefPix = [&](int x, int y) {
      if (x < 0 || x >= pixW || y < 0 || y >= pixH) return;
      const int textRow = y / 4;
      const int textCol = x / 2;
      if (inScoreBox(textRow, textCol)) return;
      refCells[textRow * plotCols + textCol] |= kBrailleBit[y % 4][x % 2];
    };
    for (int xp = 0; xp < pixW; ++xp) {
      const double xLeft  = std::max(0.0, xp - 0.5);
      const double xRight = std::min((double)(pixW - 1), xp + 0.5);
      const int yL = sampleRefYAt(xLeft);
      const int yC = sampleRefYAt(xp);
      const int yR = sampleRefYAt(xRight);
      const int yMin = std::min({yL, yC, yR});
      const int yMax = std::max({yL, yC, yR});
      for (int y = yMin; y <= yMax; ++y) setRefPix(xp, y);
    }
    Attr a(Color(CP_REFERENCE) | A_BOLD);
    char buf[4];
    for (int r = 0; r < plotRows; ++r) {
      for (int c = 0; c < plotCols; ++c) {
        const unsigned char mask = refCells[r * plotCols + c];
        if (mask == 0) continue;
        EncodeBraille(mask, buf);
        mvaddstr(plotTop + r, plotLeft + c, buf);
      }
    }
  }

  // Render braille cells (live GA curve in cyan, on top of reference).
  {
    Attr a(Color(CP_CURVE) | A_BOLD);
    char buf[4];
    for (int r = 0; r < plotRows; ++r) {
      for (int c = 0; c < plotCols; ++c) {
        const unsigned char mask = cells[r * plotCols + c];
        if (mask == 0) continue;
        EncodeBraille(mask, buf);
        mvaddstr(plotTop + r, plotLeft + c, buf);
      }
    }
  }

  // Pinned best-score box (top-right corner, drawn AFTER the curve so
  // the border sits cleanly over any cells the clip above kept empty).
  if (showScoreBox) {
    // Erase the cells the box occupies (the clip kept the curve out,
    // but axis tick rows may have written here on earlier frames).
    {
      Attr a(Color(CP_BORDER));
      for (int r = boxTop; r <= boxBottom; ++r) {
        for (int c = boxLeft; c <= boxRight; ++c) {
          mvaddstr(r, c, " ");
        }
      }
    }
    DrawBox(boxTop, boxLeft, kBoxHeight, kBoxWidth);
    const std::string boxTitle =
        bestScore.title.empty() ? "Best score" : bestScore.title;
    {
      Attr a(Color(CP_TITLE) | A_BOLD);
      const std::string titleBar = "─ " + boxTitle + " ";
      mvaddstr(boxTop, boxLeft + 2, titleBar.c_str());
    }
    char vbuf[40];
    // Custom-labelled values (mix mode's Fst) sit in [0,1] and read
    // better as fixed-point; the GA's SCval is best as %.4g.
    if (bestScore.title.empty()) {
      std::snprintf(vbuf, sizeof(vbuf), "%.4g", bestScore.score);
    } else {
      std::snprintf(vbuf, sizeof(vbuf), "%.4f", bestScore.score);
    }
    {
      Attr a(Color(CP_CURVE) | A_BOLD);
      mvaddstr(boxTop + 1, boxLeft + 2, vbuf);
    }
    char rbuf[40];
    if (bestScore.subtitle.empty()) {
      std::snprintf(rbuf, sizeof(rbuf), "round %d / %d",
                    bestScore.completedRounds, bestScore.totalRounds);
    } else {
      std::snprintf(rbuf, sizeof(rbuf), "%s", bestScore.subtitle.c_str());
    }
    {
      Attr a(Color(CP_STATUS) | A_DIM);
      mvaddstr(boxTop + 2, boxLeft + 2, rbuf);
    }
  }

  // X-axis at the bottom of the plot area.
  const int xAxisRow = plotBottom + 1;
  {
    Attr a(Color(CP_AXIS));
    mvaddstr(xAxisRow, plotLeft - 1, "└");
    for (int c = 0; c < plotCols; ++c) {
      mvaddstr(xAxisRow, plotLeft + c, "─");
    }
  }
  // X labels.
  {
    Attr a(Color(CP_AXIS_LABEL) | A_DIM);
    char lhs[16], rhs[16], mid[16];
    std::snprintf(lhs, sizeof(lhs), "gen 1");
    std::snprintf(rhs, sizeof(rhs), "gen %d", N);
    std::snprintf(mid, sizeof(mid), "gen %d", N / 2);
    mvaddstr(xAxisRow + 1, plotLeft, lhs);
    const int midCol = plotLeft + plotCols / 2 - (int)std::strlen(mid) / 2;
    mvaddstr(xAxisRow + 1, midCol, mid);
    const int rcol = std::max<int>(plotLeft, plotLeft + plotCols
                                              - (int)std::strlen(rhs));
    mvaddstr(xAxisRow + 1, rcol, rhs);
  }
}

// ---- Frame ---------------------------------------------------------------

void Render(const ProgressSnapshot& s) {
  std::lock_guard<std::mutex> lk(g_screenMu);
  erase();
  int rows = 0, cols = 0;
  getmaxyx(stdscr, rows, cols);
  if (rows < 12 || cols < 60) {
    Attr a(Color(CP_STATUS));
    mvaddstr(0, 0, "Terminal too small — resize to at least 60x12.");
    refresh();
    return;
  }

  // Title.
  DrawCenteredTitle(0, cols, "GONE2 — Genetic Optimization for Ne Estimation");

  // Global progress.
  {
    Attr a(Color(CP_STATUS));
    mvaddstr(2, kLeftMargin, "Overall");
  }
  const int gbarCol = kLeftMargin + 9;
  const int gbarRight = cols - 10;
  const int gbarWidth = std::max(10, gbarRight - gbarCol);
  DrawProgressBar(2, gbarCol, gbarWidth, s.globalProgress);
  DrawPctRight(2, cols - 8, s.globalProgress);

  // Phase list.
  const int phaseTop = 4;
  DrawPhases(phaseTop, kLeftMargin, cols, s.phases);

  // Status detail.
  const int statusRow = phaseTop + (int)s.phases.size() + 1;
  if (statusRow < rows - 2 && !s.statusDetail.empty()) {
    Attr a(Color(CP_STATUS));
    mvaddstr(statusRow, kLeftMargin, "› ");
    mvaddnstr(statusRow, kLeftMargin + 2, s.statusDetail.c_str(),
              std::min<int>(s.statusDetail.size(), cols - kLeftMargin - 4));
  }

  // Chart fills the remaining vertical space.
  const int chartTop = statusRow + 2;
  const int chartHeight = rows - chartTop - 2;
  if (chartHeight >= 8) {
    DrawNeChart(chartTop, kLeftMargin, chartHeight, cols - 2 * kLeftMargin,
                s.ne, s.bestScore);
  }

  // Hint.
  {
    Attr a(Color(CP_HINT) | A_DIM);
    const char* hint = g_finished.load(std::memory_order_relaxed)
                           ? "press Ctrl-C to exit"
                           : "press Ctrl-C to abort";
    mvaddstr(rows - 1, std::max<int>(0, cols - (int)std::strlen(hint) - 1),
             hint);
  }
  refresh();
}

void RendererLoop() {
  // Don't let SIGINT land on this thread — we want the main thread
  // to receive it so WaitForSigint can unblock cleanly.
  sigset_t set;
  sigemptyset(&set);
  sigaddset(&set, SIGINT);
  pthread_sigmask(SIG_BLOCK, &set, nullptr);

  while (!g_stop.load(std::memory_order_relaxed)) {
    ProgressSnapshot snap = g_progress->Snapshot();
    Render(snap);
    // Drain any input (ignored) so KEY_RESIZE doesn't pile up.
    int ch;
    while ((ch = getch()) != ERR) { (void)ch; }
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
  }
  if (g_progress) Render(g_progress->Snapshot());
}

extern "C" void HandleSigint(int) {
  g_sigintReceived = 1;
}

}  // namespace

extern "C" void AtExitCleanupTui() {
  // Idempotent: a normal program path that called StopNcursesTui()
  // explicitly will have already set g_running=false, so this is a
  // no-op then. On any exit() path that bypasses normal cleanup
  // (HandleInput error, ReadVcf/ReadMap failures, etc.) this runs
  // during the C++ atexit chain and restores the terminal before
  // the runtime destroys our joinable std::thread (which would
  // otherwise call std::terminate).
  StopNcursesTui();
}

void StartNcursesTui(ProgressStatus* progress) {
  if (g_running.exchange(true)) return;
  g_progress = progress;
  g_stop.store(false);

  // Capture cerr writes for the duration of ncurses mode. Restored
  // in StopNcursesTui after endwin().
  if (g_cerrOrigBuf == nullptr) {
    g_cerrCapture = new std::stringstream;
    g_cerrOrigBuf = std::cerr.rdbuf();
    std::cerr.rdbuf(g_cerrCapture->rdbuf());
  }

  // UTF-8 locale is needed for braille and box-drawing glyphs.
  std::setlocale(LC_ALL, "");
  initscr();
  cbreak();
  noecho();
  curs_set(0);
  nodelay(stdscr, TRUE);
  keypad(stdscr, TRUE);
  InitColors();

  // Make sure exit() / abort() restore the terminal. Registered
  // once per process; std::atexit guarantees the handlers run in
  // reverse registration order before static destruction.
  static bool atexitRegistered = false;
  if (!atexitRegistered) {
    std::atexit(AtExitCleanupTui);
    atexitRegistered = true;
  }

  g_thread = std::thread(RendererLoop);
}

void WaitForSigint() {
  if (!g_running.load()) return;
  g_finished.store(true, std::memory_order_relaxed);
  g_sigintReceived = 0;

  struct sigaction sa{};
  sa.sa_handler = HandleSigint;
  sigemptyset(&sa.sa_mask);
  // No SA_RESTART: we want pause() to return on the signal.
  sa.sa_flags = 0;
  struct sigaction old{};
  sigaction(SIGINT, &sa, &old);

  // Unblock SIGINT on this (main) thread in case something earlier
  // masked it.
  sigset_t set;
  sigemptyset(&set);
  sigaddset(&set, SIGINT);
  pthread_sigmask(SIG_UNBLOCK, &set, nullptr);

  while (g_sigintReceived == 0) {
    pause();
  }

  // Restore default disposition: a second Ctrl-C during cleanup will
  // terminate the process if anything hangs.
  sigaction(SIGINT, &old, nullptr);
}

void StopNcursesTui() {
  if (!g_running.exchange(false)) return;
  g_stop.store(true);
  if (g_thread.joinable()) g_thread.join();

  {
    std::lock_guard<std::mutex> lk(g_screenMu);
    endwin();
  }
  g_progress = nullptr;

  // Restore cerr and replay any error messages captured while
  // ncurses was up. After endwin() the terminal is back to its
  // pre-ncurses state, so the message lands on a clean line.
  if (g_cerrOrigBuf != nullptr) {
    std::cerr.rdbuf(g_cerrOrigBuf);
    g_cerrOrigBuf = nullptr;
    if (g_cerrCapture != nullptr) {
      const std::string captured = g_cerrCapture->str();
      if (!captured.empty()) {
        std::cerr << captured;
        std::cerr.flush();
      }
      delete g_cerrCapture;
      g_cerrCapture = nullptr;
    }
  }
}

#endif  // GONE_NCURSES_TUI

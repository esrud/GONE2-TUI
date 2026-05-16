#include "tui.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include <sys/ioctl.h>
#include <unistd.h>

namespace {

std::string Percent(double value) {
  std::ostringstream out;
  out << std::fixed << std::setprecision(2) << value * 100.0 << "%";
  return out.str();
}

std::string Bar(double progress, int width) {
  const int clampedWidth = std::max(1, width);
  const int filled = std::max(
      0, std::min(clampedWidth,
                  static_cast<int>(std::round(progress * clampedWidth))));
  return "[" + std::string(filled, '#') +
         std::string(clampedWidth - filled, '.') + "]";
}

std::string PhaseMarker(ProgressPhaseState state) {
  if (state == ProgressPhaseState::Done) return "done";
  if (state == ProgressPhaseState::Current) return "run ";
  return "wait";
}

std::vector<std::string> RenderNeChart(const NeEstimateView& ne, int width,
                                       int height) {
  std::vector<std::string> lines;
  if (!ne.hasValues || ne.values.empty() || width < 20 || height < 4) {
    lines.push_back("Waiting for Ne estimates...");
    return lines;
  }
  const double minValue =
      std::max(1.0, *std::min_element(ne.values.begin(), ne.values.end()));
  const double maxValue =
      std::max(minValue, *std::max_element(ne.values.begin(), ne.values.end()));
  const double logMin = std::log10(minValue);
  const double logMax = std::log10(maxValue);
  const double span = std::max(0.000001, logMax - logMin);
  std::vector<std::string> canvas(height, std::string(width, ' '));
  const int plotWidth = width - 10;
  for (int x = 0; x < plotWidth; ++x) {
    const std::size_t idx = static_cast<std::size_t>(std::round(
        (static_cast<double>(x) / std::max(1, plotWidth - 1)) *
        (ne.values.size() - 1)));
    const double value = std::max(1.0, ne.values[idx]);
    const double normalized = (std::log10(value) - logMin) / span;
    const int y = std::max(
        0, std::min(height - 1,
                    height - 1 -
                        static_cast<int>(std::round(normalized *
                                                    (height - 1)))));
    canvas[y][x + 9] = '*';
  }
  std::ostringstream topLabel;
  topLabel << std::fixed << std::setprecision(0) << maxValue;
  canvas.front().replace(0, std::min<std::size_t>(8, topLabel.str().size()),
                         topLabel.str().substr(0, 8));
  std::ostringstream bottomLabel;
  bottomLabel << std::fixed << std::setprecision(0) << minValue;
  canvas.back().replace(0, std::min<std::size_t>(8, bottomLabel.str().size()),
                        bottomLabel.str().substr(0, 8));
  return canvas;
}

void TerminalSize(int* width, int* height) {
  winsize size {};
  if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &size) == 0 &&
      size.ws_col > 0 && size.ws_row > 0) {
    *width = size.ws_col;
    *height = size.ws_row;
  }
}

}  // namespace

bool StdoutIsTerminal() {
  return isatty(STDOUT_FILENO) == 1;
}

std::string RenderPlainProgress(const ProgressSnapshot& snapshot) {
  std::ostringstream out;
  out << snapshot.currentPhaseName << " " << Percent(snapshot.phaseProgress)
      << " global " << Percent(snapshot.globalProgress);
  if (!snapshot.statusDetail.empty()) {
    out << " - " << snapshot.statusDetail;
  }
  out << "\n";
  return out.str();
}

std::string RenderAnsiTui(const ProgressSnapshot& snapshot, int width,
                          int height) {
  const int safeWidth = std::max(60, width);
  const int safeHeight = std::max(18, height);
  const int sideWidth = 30;
  const int mainWidth = safeWidth - sideWidth - 3;
  std::ostringstream out;
  out << "\033[?25l\033[H";
  out << "GONE2 progress" << std::string(std::max(1, safeWidth - 14), ' ')
      << "\n";
  out << "Global " << Bar(snapshot.globalProgress, mainWidth - 12) << " "
      << Percent(snapshot.globalProgress) << "\n\n";
  const std::size_t maxPhaseRows =
      static_cast<std::size_t>(std::max(1, safeHeight - 6));
  for (std::size_t row = 0; row < maxPhaseRows; ++row) {
    if (row < snapshot.phases.size()) {
      const ProgressPhaseView& phase = snapshot.phases[row];
      out << "[" << PhaseMarker(phase.state) << "] ";
      out << phase.name.substr(0, sideWidth - 8);
    } else {
      out << std::string(sideWidth, ' ');
    }
    out << " | ";
    if (row == 0) {
      out << snapshot.currentPhaseName;
    } else if (row == 1) {
      out << "Phase  " << Bar(snapshot.phaseProgress, mainWidth - 15) << " "
          << Percent(snapshot.phaseProgress);
    } else if (row == 2 && !snapshot.statusDetail.empty()) {
      out << snapshot.statusDetail.substr(0, mainWidth);
    } else if (row >= 4) {
      const std::vector<std::string> chart =
          RenderNeChart(snapshot.ne, mainWidth, safeHeight - 9);
      const std::size_t chartRow = row - 4;
      if (chartRow == 0 && snapshot.ne.hasValues) {
        out << snapshot.ne.unitLabel;
      } else if (chartRow > 0 && chartRow - 1 < chart.size()) {
        out << chart[chartRow - 1];
      }
    }
    out << "\n";
  }
  return out.str();
}

void RenderProgressToTerminal(const ProgressSnapshot& snapshot,
                              bool terminalUiEnabled,
                              std::ostream& out) {
  if (!terminalUiEnabled) return;
  if (!StdoutIsTerminal()) {
    out << RenderPlainProgress(snapshot);
    return;
  }
  int width = 100;
  int height = 30;
  TerminalSize(&width, &height);
  out << RenderAnsiTui(snapshot, width, height);
  out.flush();
}

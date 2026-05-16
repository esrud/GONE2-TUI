#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <string>
#include <vector>

enum class ProgressPhaseState {
  Pending,
  Current,
  Done
};

struct ProgressPhaseView {
  std::string name;
  ProgressPhaseState state = ProgressPhaseState::Pending;
  double progress = 0.0;
};

struct NeEstimateView {
  bool hasValues = false;
  std::vector<double> values;
  std::string unitLabel;
  // Optional ground-truth Ne by generation (1-indexed, gen 1 = most
  // recent). Loaded from -f <path>; drawn over the GA estimate in
  // red by the ncurses TUI.
  std::vector<double> referenceValues;
};

// Live "best score" indicator pinned to the top-right of the Ne chart.
//
// Two display flavours, chosen by which fields are populated:
//   - GA mode (non-mix): completedRounds/totalRounds > 0 → the TUI
//     renders "Best score / <SCval> / round N of M". This is the
//     classic GA convergence indicator.
//   - Mix mode: title/subtitle are non-empty → the TUI uses those
//     verbatim and the round counter is hidden. Lets sample.cpp put
//     the refining Fst estimate in the same corner box without GA
//     rounds existing.
struct BestScoreView {
  bool hasValue = false;
  double score = 0.0;
  int completedRounds = 0;
  int totalRounds = 0;
  std::string title;     // empty → "Best score"
  std::string subtitle;  // empty → "round N / M"
};

struct ProgressSnapshot {
  double totalTasks = 1.0;
  double currentTask = 0.0;
  double totalSteps = 1.0;
  double currentStep = 0.0;
  double phaseProgress = 0.0;
  double globalProgress = 0.0;
  std::size_t currentPhaseIndex = 0;
  std::string currentPhaseName;
  std::string statusDetail;
  std::vector<ProgressPhaseView> phases;
  NeEstimateView ne;
  BestScoreView bestScore;
};

class ProgressStatus {
 private:
  double totalTasks = 1.0;
  double currentTask = 0.0;
  double totalSteps = 1.0;
  double currentStep = 0.0;
  double taskProgress = 0.0;
  double globalProgress = 0.0;
  std::string taskName;
  std::string progressFname;
  std::string statusDetail;
  std::vector<ProgressPhaseView> phases;
  NeEstimateView ne;
  BestScoreView bestScore;
  bool terminalUiEnabled = false;
  mutable std::mutex mu;

  void EnsurePhaseList();
  void RefreshPhaseStates();
  void WriteLegacyProgressFile() const;
  void RenderTerminalProgress() const;

 public:
  void ConfigurePhases(const std::vector<std::string>& phaseNames);
  void SetTerminalUiEnabled(bool enabled);
  void InitTotalTasks(float nTasks, const char* fname);
  void SetCurrentTask(float cTask, const char* tName);
  void InitCurrentTask(float nSteps);
  void SetTaskProgress(float cStep);
  void SetStatusDetail(const std::string& detail);
  void SetNeSnapshot(const std::vector<double>& values,
                     const std::string& unitLabel);
  // Replace the optional reference Ne curve overlaid in red by the
  // ncurses TUI. The plain (non-ncurses) build silently ignores it.
  void SetReferenceNe(const std::vector<double>& values);
  // Report the best SCval seen so far. `completed` is the count of
  // GA rounds that have produced a score; `total` is the configured
  // GONE_ROUNDS so the TUI can render "round N of M".
  void SetBestScore(double score, int completed, int total);
  // Mix-mode variant: the score-box rectangle gets `title` + a value
  // line + `subtitle`, no round counter. Used to expose the refining
  // partial-Fst estimate in -x mode where there's no GA.
  void SetBestScoreLabeled(double score, const std::string& title,
                           const std::string& subtitle);
  ProgressSnapshot Snapshot() const;
  void SaveProgress();
  void PrintProgress();
};

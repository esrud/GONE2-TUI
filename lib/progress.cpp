#include "./progress.hpp"

#include <algorithm>

#ifdef GONE_TUI
#include "tui.hpp"
#endif

namespace {

double Clamp01(double value) {
  return std::max(0.0, std::min(1.0, value));
}

double PositiveOrOne(double value) {
  return value > 0.0 ? value : 1.0;
}

}  // namespace

void ProgressStatus::ConfigurePhases(
    const std::vector<std::string>& phaseNames) {
  std::lock_guard<std::mutex> lk(mu);
  phases.clear();
  phases.reserve(phaseNames.size());
  for (const std::string& name : phaseNames) {
    ProgressPhaseView phase;
    phase.name = name;
    phases.push_back(phase);
  }
  RefreshPhaseStates();
}

void ProgressStatus::SetTerminalUiEnabled(bool enabled) {
  std::lock_guard<std::mutex> lk(mu);
  terminalUiEnabled = enabled;
}

void ProgressStatus::InitTotalTasks(float nTasks, const char* fname) {
  std::lock_guard<std::mutex> lk(mu);
  totalTasks = PositiveOrOne(static_cast<double>(nTasks));
  progressFname = fname == nullptr ? "" : fname;
  EnsurePhaseList();
  RefreshPhaseStates();
}

void ProgressStatus::SetCurrentTask(float cTask, const char* tName) {
  std::unique_lock<std::mutex> lk(mu);
  if (cTask > totalTasks) {
    lk.unlock();
    std::cerr << "ERROR: Setting progress larger than max value\n";
    return;
  }
  currentTask = std::max(0.0, static_cast<double>(cTask));
  taskName = tName == nullptr ? "" : tName;
  taskProgress = 0.0;
  currentStep = 0.0;
  globalProgress = Clamp01(currentTask / totalTasks);
  EnsurePhaseList();
  RefreshPhaseStates();
}

void ProgressStatus::InitCurrentTask(float nSteps) {
  std::lock_guard<std::mutex> lk(mu);
  totalSteps = PositiveOrOne(static_cast<double>(nSteps));
  currentStep = 0.0;
  taskProgress = 0.0;
  globalProgress = Clamp01(currentTask / totalTasks);
  RefreshPhaseStates();
}

void ProgressStatus::SetTaskProgress(float cStep) {
  std::unique_lock<std::mutex> lk(mu);
  if (cStep > totalSteps) {
    lk.unlock();
    std::cerr << "ERROR: Setting progress bigger than max value\n";
    return;
  }
  currentStep = std::max(0.0, static_cast<double>(cStep));
  taskProgress = Clamp01(currentStep / totalSteps);
  globalProgress = Clamp01((currentTask + taskProgress) / totalTasks);
  RefreshPhaseStates();
  lk.unlock();
  SaveProgress();
}

void ProgressStatus::SetStatusDetail(const std::string& detail) {
  std::lock_guard<std::mutex> lk(mu);
  statusDetail = detail;
}

void ProgressStatus::SetNeSnapshot(const std::vector<double>& values,
                                   const std::string& unitLabel) {
  std::lock_guard<std::mutex> lk(mu);
  ne.hasValues = !values.empty();
  ne.values = values;
  ne.unitLabel = unitLabel;
}

void ProgressStatus::SetReferenceNe(const std::vector<double>& values) {
  std::lock_guard<std::mutex> lk(mu);
  ne.referenceValues = values;
}

void ProgressStatus::SetBestScore(double score, int completed, int total) {
  std::lock_guard<std::mutex> lk(mu);
  // Only overwrite on strict improvement so the displayed value never
  // moves backwards — the GA produces a fresh score per round and
  // they can be either better or worse than the running minimum.
  if (!bestScore.hasValue || score < bestScore.score) {
    bestScore.score = score;
  }
  bestScore.hasValue = true;
  bestScore.completedRounds = completed;
  bestScore.totalRounds = total;
  bestScore.title.clear();
  bestScore.subtitle.clear();
}

void ProgressStatus::ResetBestScore() {
  std::lock_guard<std::mutex> lk(mu);
  bestScore = BestScoreView{};
}

void ProgressStatus::SetBestScoreLabeled(double score,
                                         const std::string& title,
                                         const std::string& subtitle) {
  std::lock_guard<std::mutex> lk(mu);
  // Mix-mode value (typically a refining Fst estimate). Latest wins —
  // the metric isn't strictly monotone, and the user wants to see it
  // settle rather than stick at an early outlier.
  bestScore.score = score;
  bestScore.hasValue = true;
  bestScore.title = title;
  bestScore.subtitle = subtitle;
  bestScore.completedRounds = 0;
  bestScore.totalRounds = 0;
}

ProgressSnapshot ProgressStatus::Snapshot() const {
  std::lock_guard<std::mutex> lk(mu);
  ProgressSnapshot snapshot;
  snapshot.totalTasks = totalTasks;
  snapshot.currentTask = currentTask;
  snapshot.totalSteps = totalSteps;
  snapshot.currentStep = currentStep;
  snapshot.phaseProgress = taskProgress;
  snapshot.globalProgress = globalProgress;
  snapshot.currentPhaseIndex =
      static_cast<std::size_t>(std::max(0.0, currentTask));
  snapshot.currentPhaseName = taskName;
  snapshot.statusDetail = statusDetail;
  snapshot.phases = phases;
  snapshot.ne = ne;
  snapshot.bestScore = bestScore;
  return snapshot;
}

void ProgressStatus::SaveProgress() {
  WriteLegacyProgressFile();
  RenderTerminalProgress();
}

void ProgressStatus::PrintProgress() {
  std::ifstream progressFile(progressFname, std::ifstream::in);
  if (!progressFile.good()) return;
  std::cout << progressFile.rdbuf();
}

void ProgressStatus::EnsurePhaseList() {
  if (!phases.empty()) return;
  const int count = static_cast<int>(totalTasks);
  phases.reserve(count);
  for (int i = 0; i < count; ++i) {
    ProgressPhaseView phase;
    phase.name = "Task " + std::to_string(i + 1);
    phases.push_back(phase);
  }
}

void ProgressStatus::RefreshPhaseStates() {
  const std::size_t current =
      static_cast<std::size_t>(std::max(0.0, currentTask));
  for (std::size_t i = 0; i < phases.size(); ++i) {
    phases[i].progress = 0.0;
    if (i < current) {
      phases[i].state = ProgressPhaseState::Done;
      phases[i].progress = 1.0;
    } else if (i == current) {
      phases[i].state = ProgressPhaseState::Current;
      phases[i].progress = taskProgress;
    } else {
      phases[i].state = ProgressPhaseState::Pending;
    }
  }
}

void ProgressStatus::WriteLegacyProgressFile() const {
  if (progressFname.empty()) return;
  std::ofstream progressFile(progressFname, std::ios::out);
  if (!progressFile.good()) return;
  progressFile << "    GLOBAL PROGRESS: " << std::fixed
               << std::setprecision(2) << globalProgress * 100.0 << "%\n";
  for (const ProgressPhaseView& phase : phases) {
    if (phase.state == ProgressPhaseState::Done) {
      progressFile << "    " << phase.name << " ... DONE\n";
    }
  }
  progressFile << "    " << taskName << " ... "
               << taskProgress * 100.0 << "%\n";
}

void ProgressStatus::RenderTerminalProgress() const {
#ifdef GONE_TUI
  RenderProgressToTerminal(Snapshot(), terminalUiEnabled, std::cout);
#endif
}

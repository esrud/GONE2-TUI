#pragma once

#include <iosfwd>
#include <string>

#include "progress.hpp"

bool StdoutIsTerminal();
std::string RenderPlainProgress(const ProgressSnapshot& snapshot);
std::string RenderAnsiTui(const ProgressSnapshot& snapshot, int width,
                          int height);
void RenderProgressToTerminal(const ProgressSnapshot& snapshot,
                              bool terminalUiEnabled,
                              std::ostream& out);

#pragma once

#include <string>
#include <vector>

// Parse a ground-truth Ne history file. Two accepted formats:
//   1. simu params (trailing "indiv_diploides generaciones" block of
//      `Ne nGen` rows, oldest first); returns one Ne value per
//      generation, gen 1 = most recent.
//   2. per-generation "Gen Ne" table (data_*.true). Gen indices may
//      skip; gaps are constant-extended from the previous Ne.
// Returns an empty vector if the file can't be parsed.
std::vector<double> LoadReferenceNe(const std::string& path);

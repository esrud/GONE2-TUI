#include "ne_ref.hpp"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

// Strip leading whitespace and trailing CR/whitespace in place.
void Trim(std::string& s) {
  while (!s.empty() && (s.front() == ' ' || s.front() == '\t')) s.erase(0, 1);
  while (!s.empty() && (s.back() == ' ' || s.back() == '\t' ||
                        s.back() == '\r' || s.back() == '\n')) s.pop_back();
}

// True if every char in s is a digit (s non-empty).
bool IsAllDigits(const std::string& s) {
  if (s.empty()) return false;
  for (char c : s) if (c < '0' || c > '9') return false;
  return true;
}

}  // namespace

std::vector<double> LoadReferenceNe(const std::string& path) {
  std::ifstream in(path);
  if (!in.good()) {
    std::cerr << "Could not open Ne reference file: " << path << "\n";
    return {};
  }

  // Block format: ignore everything until the marker line, then read
  // pairs of `Ne nGen` (oldest → most recent), expand in reverse.
  // Per-gen format: read all 2-column numeric rows, sort by Gen.
  std::vector<std::pair<long long, long long>> blocks;
  std::vector<std::pair<long long, double>> perGen;
  bool inBlock = false;
  std::string line;
  while (std::getline(in, line)) {
    Trim(line);
    if (line.empty()) continue;
    if (!line.empty() && line[0] == '#') {
      // Comment. Watch for the simu-params block marker.
      if (line.find("EL PRIMER BLOQUE") != std::string::npos) {
        inBlock = true;
      }
      continue;
    }
    // Read up to two tokens.
    std::istringstream iss(line);
    std::string t1, t2;
    if (!(iss >> t1 >> t2)) continue;
    // Per-gen format header: "Gen Ne" — non-numeric tokens.
    if (!std::isdigit(static_cast<unsigned char>(t1[0])) &&
        t1[0] != '-' && t1[0] != '+') {
      continue;
    }
    if (inBlock && IsAllDigits(t1) && IsAllDigits(t2)) {
      blocks.emplace_back(std::stoll(t1), std::stoll(t2));
    } else {
      // Per-gen: first column must be a non-negative integer gen, second
      // a positive number for Ne.
      try {
        const long long gen = std::stoll(t1);
        const double ne     = std::stod(t2);
        if (gen >= 0 && ne > 0) perGen.emplace_back(gen, ne);
      } catch (...) {
        // ignore
      }
    }
  }

  std::vector<double> out;
  if (!blocks.empty()) {
    // Block format. The last entry is most recent.
    for (auto it = blocks.rbegin(); it != blocks.rend(); ++it) {
      const double ne = static_cast<double>(it->first);
      const long long n_gen = it->second;
      for (long long i = 0; i < n_gen; ++i) out.push_back(ne);
    }
    return out;
  }
  if (!perGen.empty()) {
    std::sort(perGen.begin(), perGen.end());
    long long lastGen = -1;
    double lastNe = perGen[0].second;
    for (const auto& [gen, ne] : perGen) {
      while (lastGen + 1 < gen) {
        ++lastGen;
        out.push_back(lastNe);
      }
      lastGen = gen;
      lastNe = ne;
      out.push_back(ne);
    }
    return out;
  }
  std::cerr << "Could not parse Ne reference from: " << path
            << " (expected simu params or `Gen Ne` table)\n";
  return {};
}

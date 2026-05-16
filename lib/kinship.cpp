#include "kinship.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <omp.h>

#include "progress.hpp"

namespace {

// Pairwise count accumulators for the KING-robust φ at a single
// (i, j) pair, over all loci where both individuals are called.
struct PairCounts {
  long long NAaAa = 0;  // both i and j heterozygous
  long long NAAaa = 0;  // opposite homozygotes (one AA, one aa)
  long long NAaI  = 0;  // i heterozygous (and j called)
  long long NAaJ  = 0;  // j heterozygous (and i called)
};

// Compute φ from the counts. Returns 0 if both individuals have no
// heterozygous calls in common-called loci — the formula's denominator
// goes to zero and we can't make a useful estimate.
double PhiFromCounts(const PairCounts& c) {
  const double denom = static_cast<double>(c.NAaI + c.NAaJ);
  if (denom <= 0) return 0.0;
  return (static_cast<double>(c.NAaAa) -
          2.0 * static_cast<double>(c.NAAaa)) / denom;
}

}  // namespace

// Compute KING-robust φ for every (i, j) pair, returning a sparse
// upper-triangle keyed by i. Includes all pairs whose φ exceeds
// `minPhiKept` so the caller can both filter and report on tiers
// without recomputing.
static std::vector<std::vector<std::pair<int, double>>>
ComputeAllPhiPairs(const PopulationInfo* popInfo, double minPhiKept) {
  const int N = popInfo->numIndi;
  const int L = popInfo->numLoci;
  std::vector<std::vector<std::pair<int, double>>> out(N);
#pragma omp parallel for schedule(dynamic, 4)
  for (int i = 0; i < N; ++i) {
    const uint8_t* gi = popInfo->indi[i];
    std::vector<std::pair<int, double>> local;
    for (int j = i + 1; j < N; ++j) {
      const uint8_t* gj = popInfo->indi[j];
      PairCounts c;
      for (int locus = 0; locus < L; ++locus) {
        const uint8_t a = gi[locus];
        const uint8_t b = gj[locus];
        if (a >= 9 || b >= 9) continue;
        if (a == 1 && b == 1)        ++c.NAaAa;
        else if ((a == 0 && b == 2) || (a == 2 && b == 0)) ++c.NAAaa;
        if (a == 1) ++c.NAaI;
        if (b == 1) ++c.NAaJ;
      }
      const double phi = PhiFromCounts(c);
      if (phi > minPhiKept) local.emplace_back(j, phi);
    }
    out[i] = std::move(local);
  }
  return out;
}

SibshipReport ComputeSibshipReport(
    int N,
    const std::vector<std::vector<std::pair<int, double>>>& phiPairs,
    ProgressStatus* progress) {
  // KING tier cutoffs from Manichaikul et al. 2010.
  constexpr double kFirstHi  = 0.354;
  constexpr double kFirstLo  = 0.177;
  constexpr double kSecondLo = 0.0884;
  constexpr double kThirdLo  = 0.0442;
  SibshipReport r{};
  r.n = N;
  for (const auto& row : phiPairs) {
    for (const auto& kv : row) {
      const double phi = kv.second;
      if      (phi >  kFirstLo  && phi <= kFirstHi)  ++r.firstDegree;
      else if (phi >  kSecondLo && phi <= kFirstLo)  ++r.secondDegree;
      else if (phi >  kThirdLo  && phi <= kSecondLo) ++r.thirdDegree;
    }
  }
  const double pairs = (N > 1) ? static_cast<double>(N) * (N - 1) / 2.0 : 1.0;
  // Idealised Nb estimates from pair counts. Under a Wright-Fisher
  // population with Nb breeders, the probability that two random
  // sampled individuals are full sibs is ≈ 1/Nb (monogamy) or scales
  // like 2/Nb (random mating). C(n,2) / count_FS hits the right
  // order of magnitude. Half-sib version uses the same logic with a
  // ×4 fudge factor for the second-degree relationship. Both should
  // be treated as rough — the formulas assume one offspring per
  // parent pair and ignore population structure.
  r.nbFromFirst  =
      r.firstDegree  > 0 ? pairs / static_cast<double>(r.firstDegree)  : 0.0;
  r.nbFromSecond =
      r.secondDegree > 0
          ? 4.0 * pairs / static_cast<double>(r.secondDegree)
          : 0.0;

  if (progress != nullptr) {
    std::string msg =
        "Sibship: 1st-deg=" + std::to_string(r.firstDegree) +
        " 2nd-deg=" + std::to_string(r.secondDegree) +
        " 3rd-deg=" + std::to_string(r.thirdDegree);
    if (r.nbFromFirst > 0) {
      msg += "; Nb≈" + std::to_string(static_cast<long>(r.nbFromFirst)) +
             " (from 1st-deg)";
    } else if (r.nbFromSecond > 0) {
      msg += "; Nb≈" + std::to_string(static_cast<long>(r.nbFromSecond)) +
             " (from 2nd-deg)";
    }
    progress->SetStatusDetail(msg);
    std::cerr << " " << msg << "\n";
  }
  return r;
}

int FilterRelatedIndividuals(PopulationInfo* popInfo,
                             double phiThreshold,
                             ProgressStatus* progress) {
  if (phiThreshold <= 0.0) return 0;
  const int N = popInfo->numIndi;
  if (N < 2) return 0;

  if (progress != nullptr) {
    progress->SetStatusDetail(
        "Computing pairwise kinship (KING-robust) over " +
        std::to_string(N) + " individuals");
  }

  // Compute the full φ matrix once, then re-use it for both the
  // sibship report and the related-pair filter. Keep everything
  // above the 3rd-degree boundary (0.0442) — the report bins on
  // higher cutoffs and the filter uses phiThreshold ≥ 0.0442 in
  // practice, so dropping anything below saves memory and lookups.
  constexpr double kKeepFloor = 0.0442;
  const double floorPhi = std::min(kKeepFloor, phiThreshold);
  auto phiPairs = ComputeAllPhiPairs(popInfo, floorPhi);

  // Sibship-based Nb report (uses the full kept set of tiers, not
  // just edges above the filter threshold).
  ComputeSibshipReport(N, phiPairs, progress);

  // Now build the filter adjacency from pairs above phiThreshold.
  std::vector<std::vector<int>> related(N);
  for (int i = 0; i < N; ++i) {
    for (const auto& kv : phiPairs[i]) {
      if (kv.second > phiThreshold) {
        related[i].push_back(kv.first);
        related[kv.first].push_back(i);
      }
    }
  }

  // Greedy max-IS: repeatedly drop the surviving vertex with the
  // highest degree, breaking edges on its neighbours. Equivalent to
  // the standard "minimum vertex cover" heuristic used by KING and
  // SNPRelate when they prune for unrelated sets.
  std::vector<bool> dropped(N, false);
  while (true) {
    int worst = -1;
    int worstDeg = 0;
    for (int i = 0; i < N; ++i) {
      if (dropped[i]) continue;
      int deg = 0;
      for (int nb : related[i]) if (!dropped[nb]) ++deg;
      if (deg > worstDeg) { worstDeg = deg; worst = i; }
    }
    if (worst < 0) break;  // no edges left
    dropped[worst] = true;
  }

  const int nDropped =
      static_cast<int>(std::count(dropped.begin(), dropped.end(), true));
  if (nDropped == 0) {
    const std::string msg =
        "Kinship filter (-k " + std::to_string(phiThreshold) +
        "): 0 individuals dropped (no pair above threshold)";
    if (progress != nullptr) progress->SetStatusDetail(msg);
    std::cerr << " " << msg << "\n";
    return 0;
  }

  // Compact popInfo->indi[] so survivors occupy 0..numIndi-nDropped-1.
  // Genotypes are POD bytes so we can just memmove rows.
  int write = 0;
  for (int read = 0; read < N; ++read) {
    if (dropped[read]) continue;
    if (write != read) {
      std::copy(popInfo->indi[read],
                popInfo->indi[read] + popInfo->numLoci,
                popInfo->indi[write]);
    }
    ++write;
  }
  popInfo->numIndi = write;

  const std::string msg =
      "Kinship filter (-k " + std::to_string(phiThreshold) +
      "): dropped " + std::to_string(nDropped) + " of " +
      std::to_string(N) + " individuals; " +
      std::to_string(popInfo->numIndi) + " remain";
  if (progress != nullptr) progress->SetStatusDetail(msg);
  std::cerr << " " << msg << "\n";
  return nDropped;
}

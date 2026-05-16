#pragma once

#include <utility>
#include <vector>

#include "population.hpp"

class ProgressStatus;

// Drop individuals from `popInfo` whose pairwise KING-robust kinship
// coefficient exceeds `phiThreshold`. KING-robust φ is the standard
// between-family estimator from Manichaikul et al. (2010):
//
//   φ_ij = (N_Aa,Aa  −  2·N_AA,aa)  /  (N_Aa^i + N_Aa^j)
//
// where the counts are over called-in-both loci. Typical cutoffs:
//   0.354 — duplicates / MZ twins
//   0.177 — 1st-degree (parent/offspring, full sibs)
//   0.0884 — 2nd-degree (half-sibs, grandparent)
//   0.0442 — 3rd-degree (first cousins)
// Threshold 0 (or negative) is a no-op.
//
// Greedy max-independent-set on the related-pair graph: repeatedly
// drop the individual with the most surviving edges until no edges
// remain. Compacts popInfo->indi[] in place and adjusts numIndi so
// the rest of the pipeline sees the filtered sample directly.
//
// Returns the number of individuals dropped.
int FilterRelatedIndividuals(PopulationInfo* popInfo,
                             double phiThreshold,
                             ProgressStatus* progress);

// Sibship-based estimate of the effective number of breeders (Nb)
// in the parental generation. Counts close-relative pairs in each
// KING-robust φ tier and back-solves Nb under an idealised
// Wright-Fisher assumption. Reports to `progress` (and stderr in
// non-TUI builds). Caller must have already computed the pairwise
// φ matrix; we accept it as an N×N upper-triangle stored in a
// vector of vectors to avoid recomputing.
//
// `phiPairs[i]` is a list of (j, φ_ij) for j > i.
struct SibshipReport {
  int n;                 // sample size at the time of counting
  long firstDegree;      // pairs with 0.177 < φ ≤ 0.354  (FS / PO)
  long secondDegree;     // pairs with 0.0884 < φ ≤ 0.177 (HS / GP)
  long thirdDegree;      // pairs with 0.0442 < φ ≤ 0.0884 (1C)
  double nbFromFirst;    // C(n,2) / max(1, firstDegree)
  double nbFromSecond;   // 4·C(n,2) / max(1, secondDegree)
};

SibshipReport ComputeSibshipReport(
    int N, const std::vector<std::vector<std::pair<int, double>>>& phiPairs,
    ProgressStatus* progress);


#include "./bicho.hpp"

#include <cmath>

// Bin-weighting and smoothness regularisation knobs.
//
// Defaults on this experimental branch:
//   GA_SMOOTH_LAMBDA = 1e-4    smoothness regulariser strength
//   GA_SMOOTH_L1     = undef   quadratic (L2) penalty on log10(Ne)
//                              jumps. The 19-sim sweep at full
//                              strength (sweep_out/summary.tsv) shows
//                              L2 generalises best: lowest mean RMSE,
//                              smallest worst-case, most wins, and
//                              detects more bottleneck features than
//                              L1 — including BOT4 where L1 blunts
//                              the dip.
//   GA_BIN_WEIGHTED  = on      residuals weighted by per-bin SNP-
//                              pair count
//
// Build with -DGA_SMOOTH_L1=1 to use L1 (total-variation) instead.
// Build with -DGA_SMOOTH_CUTOFF=<float> for a truncated quadratic.
// See docs/superpowers/specs/2026-05-13-ga-loss-tweaks.md and
// 2026-05-13-ga-smooth-shape.md for the supporting sweep.
#ifndef GA_SMOOTH_LAMBDA
#define GA_SMOOTH_LAMBDA 0.0001
#endif
#ifndef GA_BIN_WEIGHTED
#define GA_BIN_WEIGHTED 1
#endif

inline double powah(double base, int exponent) {
  double result = 1;
  while (exponent > 0) {
    if (exponent & 1) {
      result *= base;
    }
    base *= base;
    exponent >>= 1;
  }
  return result;
}

// Core of computation:
// Calculates the predicted d²_c for a piecewise-constant Ne history
// described by `bicho` and accumulates the squared residual against
// the observed d²_c stored in `sampleInfo`.
//
// Within one call, every segment-only quantity (Ne12base, Ne122base,
// Ne102base, their powers over the segment width, and b = (1-Ne12^w)/
// (1-Ne12base)) is invariant across recombination bins. We compute
// them once into local arrays and reuse them inside the per-bin loop;
// only the c-dependent powers of (1-c)² remain per-bin work.
static double CalculaSCImpl(Bicho* bicho, GsampleInfo* sampleInfo,
                            double* d2cPred, double cutoff) {
  const int nsegmentos = bicho->nSeg;
  const double Necons = bicho->NeBl[nsegmentos - 1];

  int hastasegmento = nsegmentos;
  if (hastasegmento > 1) {
    --hastasegmento;
  }

  // Per-segment, bin-invariant precomputation.
  int    segExp[kNumLinMax];
  double segNeBl[kNumLinMax];
  double Ne12Base[kNumLinMax];
  double Ne122Base[kNumLinMax];
  double Ne102Base[kNumLinMax];
  double Ne12Pow[kNumLinMax];
  double Ne122Pow[kNumLinMax];
  double Ne102Pow[kNumLinMax];
  double bConst[kNumLinMax];   // (1 - Ne12^w) / (1 - Ne12Base)
  double swCoeff[kNumLinMax];  // (1 - Ne12^w) / (2/NeBl)  for Sw accumulation

  for (int i = 0; i < hastasegmento; ++i) {
    const int   w   = bicho->segBl[i + 1] - bicho->segBl[i];
    const double Ne = bicho->NeBl[i];
    segExp[i]    = w;
    segNeBl[i]   = Ne;
    Ne12Base[i]  = 1.0 - 2.0 / Ne;
    Ne122Base[i] = 1.0 - 2.2 / Ne;
    Ne102Base[i] = 1.0 - 0.2 / Ne;
    Ne12Pow[i]   = powah(Ne12Base[i],  w);
    Ne122Pow[i]  = powah(Ne122Base[i], w);
    Ne102Pow[i]  = powah(Ne102Base[i], w);
    bConst[i]    = (1.0 - Ne12Pow[i]) / (1.0 - Ne12Base[i]);
    swCoeff[i]   = (1.0 - Ne12Pow[i]) / (2.0 / Ne);
  }

  const double Ne12nSegMinusOne_base = 1.0 - 2.0 / Necons;
  const double oneMinus2_2overNecons = 1.0 - 2.2 / Necons;
  const double Necons_div2 = Necons / 2.0;

  const int haplotype = sampleInfo->haplotype;
  const bool mix      = sampleInfo->mix;

  double score = 0;

  // Smoothness regularisation: penalise log10(Ne) jumps between
  // adjacent segments. Three penalty shapes are available; the choice
  // is read at runtime from sampleInfo->gaConfig so a single binary
  // can sweep L1, L2, and truncated in sequence (gone-ncurses_combo).
  // Non-combo builds initialise gaConfig from the compile-time
  // -DGA_* macros and behave bit-identically to the previous code
  // path. Setting smoothLambda to 0 disables the term.
  {
    const GAConfig& cfg = sampleInfo->gaConfig;
    if (cfg.smoothLambda > 0.0) {
      double smooth = 0;
      const double truncSq = cfg.smoothCutoff * cfg.smoothCutoff;
      for (int i = 1; i < nsegmentos; ++i) {
        const double a = std::max(bicho->NeBl[i - 1], MIN_NE_SIZE);
        const double b = std::max(bicho->NeBl[i],     MIN_NE_SIZE);
        const double d = std::log10(b) - std::log10(a);
        switch (cfg.smoothKind) {
          case GASmoothKind::L1:
            smooth += std::fabs(d);
            break;
          case GASmoothKind::Truncated: {
            const double dSq = d * d;
            smooth += (dSq < truncSq) ? dSq : truncSq;
            break;
          }
          case GASmoothKind::Quadratic:
          default:
            smooth += d * d;
            break;
        }
      }
      score += cfg.smoothLambda * smooth;
    }
  }

  // Heterozygosity anchor: softly pull the ancient plateau (Necons,
  // the deepest segment's Ne) toward the mutation-drift-balance
  // target the user supplied via -m. Disabled when either field is
  // zero. Penalty is on log10 distance so the term has comparable
  // magnitude to the smoothness term at all Ne scales.
  {
    const GAConfig& cfg = sampleInfo->gaConfig;
    if (cfg.anchorLambda > 0.0 && cfg.anchorNeHaploid > 0.0) {
      const double Necons_safe = std::max(Necons, MIN_NE_SIZE);
      const double dlog = std::log10(Necons_safe) -
                          std::log10(cfg.anchorNeHaploid);
      score += cfg.anchorLambda * dlog * dlog;
    }
  }

#ifdef GA_BIN_WEIGHTED
  // Average nBin so the weighted score keeps the same magnitude as
  // the unweighted residual sum — bins with above-average pair
  // counts dominate, bins with very few pairs (large variance) are
  // down-weighted.
  double sumNBinAll = 0;
  for (int ii = 0; ii < sampleInfo->nBins; ++ii) {
    sumNBinAll += sampleInfo->nBin[ii];
  }
  const double invMeanNBin =
      (sumNBinAll > 0) ? (sampleInfo->nBins / sumNBinAll) : 1.0;
#endif

  for (int ii = 0; ii < sampleInfo->nBins; ++ii) {
    const double cv = sampleInfo->oneMinuscValSq[ii];
    double p1a = 1, p1b = 1, r1a = 1, s1 = 0;
    double Sd2 = 0, Sw = 0;
    double acuOneMinusCvalSq = 1;

    const double Nec122nSegMinusOne = cv * oneMinus2_2overNecons;

    for (int i = 0; i < hastasegmento; ++i) {
      const double Ne12ancho  = Ne12Pow[i];
      const double Ne122ancho = Ne122Pow[i];
      const double Ne102ancho = Ne102Pow[i];
      const double oneMinuscValSqAncho =
          sampleInfo->oneMinuscValSqPow[ii][segExp[i]];

      const double Ne12base  = Ne12Base[i];
      const double Ne122base = Ne122Base[i];
      const double Nec122    = cv * Ne122base;
      const double Nec102    = cv * Ne102Base[i];

      const double a = (1.0 - Ne122ancho * oneMinuscValSqAncho) / (1.0 - Nec122);
      const double b = bConst[i];

      Sd2 += s1 * p1a * b + p1b / segNeBl[i] * acuOneMinusCvalSq *
                                (Ne12base * b - Nec122 * a) /
                                (Ne12base - Nec122);
      Sw += p1a * swCoeff[i];

      s1 += r1a * acuOneMinusCvalSq / segNeBl[i] *
            (1.0 - Ne102ancho * oneMinuscValSqAncho) / (1.0 - Nec102);
      r1a *= Ne102ancho;
      p1a *= Ne12ancho;
      p1b *= Ne122ancho;
      acuOneMinusCvalSq *= oneMinuscValSqAncho;
    }

    const double aPlus = Nec122nSegMinusOne / (1.0 - Nec122nSegMinusOne);
    const double bPlus = Ne12nSegMinusOne_base * Necons_div2;

    Sd2 += s1 * p1a * Necons_div2 + p1b / Necons * acuOneMinusCvalSq *
                                       (bPlus - aPlus) /
                                       (Ne12nSegMinusOne_base - Nec122nSegMinusOne);
    Sw += p1a * Necons_div2;

    double d2c = Sd2 / Sw;

    // Successive-generations sampling correction.
    d2c *= sampleInfo->ngensamplingcorrec[ii];
    // Diploid accumulation correction × (1 + c²).
    if (haplotype != 1) {
      d2c *= sampleInfo->onePluscValSq[ii];
    }

    // Sampling correction by haplotype.
    double pred;
    switch (haplotype) {
      case 1:  // Haploids.
      case 2:  // Phased diploids.
        pred = d2c * sampleInfo->basecallcorrec * sampleInfo->cValRep[ii] +
               sampleInfo->sampleZ4;
        break;
      case 3:  // Pseudohaploids (low coverage).
        pred = d2c / 4 * sampleInfo->basecallcorrec * sampleInfo->cValRep[ii] +
               sampleInfo->sampleZ3;
        break;
      default:  // Unphased diploids.
        if (mix) {
          pred = d2c * sampleInfo->basecallcorrec;
        } else {
          pred = (d2c * sampleInfo->basecallcorrec * sampleInfo->cValRep[ii] *
                      sampleInfo->sampleX +
                  sampleInfo->cValSq[ii] / Necons * sampleInfo->sampleX +
                  sampleInfo->sampleY) /
                 sampleInfo->correccion;
        }
        break;
    }
    if (d2cPred != nullptr) {
      d2cPred[ii] = pred;
    }

    const double residual = sampleInfo->d2cObs[ii] - pred;
#ifdef GA_BIN_WEIGHTED
    const double binWeight = sampleInfo->nBin[ii] * invMeanNBin;
    score += binWeight * Square<double>(residual);
#else
    score += Square<double>(residual);
#endif
    if (d2cPred == nullptr && score > cutoff) {
      return cutoff;
    }
  }
  return score;
}

double CalculaSC(Bicho* bicho, GsampleInfo* sampleInfo, double* d2cPred) {
  return CalculaSCImpl(bicho, sampleInfo, d2cPred, MAX_DOUBLE);
}

double CalculaSCScoreOnly(Bicho* bicho, GsampleInfo* sampleInfo) {
  return CalculaSCImpl(bicho, sampleInfo, nullptr, MAX_DOUBLE);
}

double CalculaSCScoreCutoff(Bicho* bicho, GsampleInfo* sampleInfo,
                            double cutoff) {
  return CalculaSCImpl(bicho, sampleInfo, nullptr, cutoff);
}

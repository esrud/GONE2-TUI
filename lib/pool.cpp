
#include "./pool.hpp"

void SetMutationParameters(Pool* pool, double frecMut, double efectoMut,
                           double efectoMutLateral, double frecRec,
                           double frecNoRec, double frecMutLateral,
                           int maxSegmentos, bool solape) {
  pool->mutParams.frecMut          = frecMut;
  pool->mutParams.efectoMut        = efectoMut;
  pool->mutParams.efectoMutLateral = efectoMutLateral;
  pool->mutParams.frecRec          = frecRec;
  pool->mutParams.frecNoRec        = frecNoRec;
  pool->mutParams.frecMutLateral   = frecMutLateral;
  pool->mutParams.maxSegmentos     = maxSegmentos;
  pool->mutParams.solape           = solape;
}

struct ChildWorstState {
  double score;
  int index;
  bool initialized;
};

static void RefreshWorstChild(Pool* pool, ChildWorstState* worst) {
  double maxSC = -1;
  int idxMaxSC = 0;
  for (int i = 0; i < kNHijos; ++i) {
    if (pool->children[i].SCval > maxSC) {
      maxSC = pool->children[i].SCval;
      idxMaxSC = i;
    }
  }
  worst->score = maxSC;
  worst->index = idxMaxSC;
  worst->initialized = true;
}

void SetInitialPoolParameters(Pool* pool, double mincval) {
  // Maximum generation index considered by the GA. Bounded so the
  // segBl-indexed arrays never overflow.
  const int gmax = std::min(static_cast<int>(1.0 / mincval), kNumGenMax - 10);
  const int gmax2 = gmax * 26 / 40;
  const int topeposigen = gmax2 - 2 * kResolucion;

  pool->poolParams.gmax[0] = gmax;
  pool->poolParams.gmax[1] = gmax2;
  pool->poolParams.gmax[2] = gmax2 + (gmax - gmax2) / 2;

  pool->poolParams.topePosInicio[0] = std::min(topeposigen, kTopePosInicio1);
  pool->poolParams.topePosInicio[1] = std::min(topeposigen, kTopePosInicio2);
  pool->poolParams.topePosInicio[2] = std::min(topeposigen, kTopePosInicio3);
  pool->poolParams.topePosInicio[3] = std::min(topeposigen, kTopePosInicio4);

  pool->numGenerations = kNumGenerations;
  pool->numChildren    = kNumDes;

  pool->mutParams.frecMut          = 0.3;
  pool->mutParams.efectoMut        = 0.3;
  pool->mutParams.efectoMutLateral = 5.0;
  pool->mutParams.frecRec          = 0.6;
  pool->mutParams.frecNoRec        = 0.3;
  pool->mutParams.frecMutLateral   = 0.2;
  pool->mutParams.maxSegmentos     = 4;
  pool->mutParams.solape           = true;
}

void PrePopulatePool(Pool* pool, GsampleInfo* sampleInfo) {
  // Random initialisation of kNumDesInicio=3000 candidate bichos; the
  // best kTercio12=100 become the starting parent pool.
  Bicho bicho = {0};
  bicho.efval     = sampleInfo->fVal;
  bicho.nSeg      = 4;
  bicho.segBl[0]  = 0;
  bicho.segBl[3]  = pool->poolParams.gmax[2];
  bicho.segBl[4]  = pool->poolParams.gmax[0];

  double r;
  double tope, Ne, NeCoeff;
  int posi1, posi2;

  for (int i = 0; i < kNumDesInicio; ++i) {
    // Sample two block boundary positions from the same range bin
    // until they are sufficiently far apart.
    bool isValidChild = false;
    while (!isValidChild) {
      r = xprng.uniforme01();
      int rangeIdx;
      if      (r < 0.4) rangeIdx = 1;
      else if (r < 0.8) rangeIdx = 2;
      else              rangeIdx = 3;
      const int range = pool->poolParams.topePosInicio[rangeIdx] - kResolucion;
      posi1 = static_cast<int>(xprng.uniforme01() * range) + kResolucion;
      posi2 = static_cast<int>(xprng.uniforme01() * range) + kResolucion;

      if (posi1 > posi2) std::swap(posi1, posi2);
      if ((posi2 - posi1) > MIN_NE_SIZE) isValidChild = true;
    }
    bicho.segBl[1] = posi1;
    bicho.segBl[2] = posi2;
    bicho.NeBl[0]  = sampleInfo->NeMed;

    // Single amplitude `tope` for all three Ne perturbations of this
    // bicho — the first two share it, the last halves it.
    tope = 2 * xprng.uniforme01() * kTopeSalto;
    for (int j = 1; j < 3; ++j) {
      NeCoeff = 1.0 + xprng.uniforme01() * tope;
      if (xprng.uniforme01() < 0.66) NeCoeff = 1.0 / NeCoeff;
      Ne = std::max(std::min(sampleInfo->NeMed * NeCoeff, MAX_NE_SIZE),
                    MIN_NE_SIZE);
      bicho.NeBl[j] = Ne;
    }
    NeCoeff = 1.0 + xprng.uniforme01() * tope / 2.0;
    Ne = std::max(std::min(sampleInfo->NeMed * NeCoeff, MAX_NE_SIZE),
                  MIN_NE_SIZE);
    bicho.NeBl[3] = Ne;
    if (i < kTercio12) {
      bicho.SCval = CalculaSCScoreOnly(&bicho, sampleInfo);
      pool->parents[i] = bicho;
    } else {
      // Replace the worst parent if this candidate is better.
      double maxSC = 0;
      int idxMaxSC = 0;
      for (int j = 0; j < kTercio12; ++j) {
        if (pool->parents[j].SCval > maxSC) {
          idxMaxSC = j;
          maxSC = pool->parents[j].SCval;
        }
      }
      bicho.SCval = CalculaSCScoreCutoff(&bicho, sampleInfo, maxSC);
      if (bicho.SCval < pool->parents[idxMaxSC].SCval) {
        pool->parents[idxMaxSC] = bicho;
      }
    }
  }

  // Sort parents by ascending SCval. NOTE: this is the original
  // selection sort with a quirk — it does not consider parents[i]
  // when picking the minimum, so the result differs from a textbook
  // selection sort. Preserved verbatim because changing it shifts
  // the GA into a different basin and changes Ne outputs.
  for (int i = 0; i < kTercio12 - 1; ++i) {
    int idxMinSC = i;
    double minSC = MAX_DOUBLE;
    for (int j = i + 1; j < kTercio12; ++j) {
      if (pool->parents[j].SCval < minSC) {
        idxMinSC = j;
        minSC = pool->parents[j].SCval;
      }
    }
    if (idxMinSC != i) {
      bicho = pool->parents[i];
      pool->parents[i] = pool->parents[idxMinSC];
      pool->parents[idxMinSC] = bicho;
    }
  }
}

// Reject if any existing segment boundary is closer than kResolucion
// generations to the proposed crossover point posigen.
bool CheckParentMutations(Bicho* parent, const int posigen) {
  for (int j = 0; j <= parent->nSeg; ++j) {
    // Branchless abs for 32-bit int via sign mask.
    int abs_val = parent->segBl[j] - posigen;
    const int temp = abs_val >> 31;
    abs_val ^= temp;
    abs_val += temp & 1;
    if (abs_val <= kResolucion) return false;
  }
  return true;
}

// Try to fuse two adjacent blocks into one. Returns false if the
// resulting Ne ratio with a neighbouring block would exceed the
// allowed jump factor, or if the segment count falls outside bounds.
bool MutateFusion(Bicho* bicho, Pool* pool) {
  if (bicho->nSeg > 3) {
    if ((xprng.uniforme01() < pool->mutParams.frecNoRec) ||
        (bicho->nSeg > pool->mutParams.maxSegmentos)) {
      const int posiblock =
          static_cast<int>(xprng.uniforme01() * (bicho->nSeg - 2));

      const double ancho1 =
          bicho->segBl[posiblock + 1] - bicho->segBl[posiblock];
      const double ancho2 =
          bicho->segBl[posiblock + 2] - bicho->segBl[posiblock + 1];

      const double ancho =
          (ancho1 + ancho2) /
          (ancho1 / bicho->NeBl[posiblock] +
           ancho2 / bicho->NeBl[posiblock + 1]);
      if (posiblock > 0 &&
          ancho / bicho->NeBl[posiblock - 1] < kInvTopeSalto2) {
        return false;
      }
      if (posiblock < bicho->nSeg - 2 &&
          ancho / bicho->NeBl[posiblock + 2] > kTopeSalto2) {
        return false;
      }

      bicho->NeBl[posiblock] = ancho;
      for (int j = posiblock + 1; j <= bicho->nSeg; ++j) {
        bicho->NeBl[j]  = bicho->NeBl[j + 1];
        bicho->segBl[j] = bicho->segBl[j + 1];
      }
      --bicho->nSeg;
    }
  }
  if (bicho->nSeg > pool->mutParams.maxSegmentos || bicho->nSeg < 3) {
    return false;
  }
  return true;
}

void Reproduce(Pool* pool, GsampleInfo* sampleInfo, int childIdx,
               ChildWorstState* worstChild) {
  // Build a single new bicho via either recombination of two parents
  // (50% inter-parent + 50% intra-parent) or straight copy of parent1.
  Bicho child = {0};
  Bicho* parent1 = nullptr;
  Bicho* parent2 = nullptr;
  bool isValidChild = false;
#if FLAG_POR_BLOQUES
  int posiblock;
#endif
  int posigen = 0;

  while (!isValidChild) {
    int ind1 = static_cast<int>(xprng.uniforme01() * kTercio1);
    int ind2 = static_cast<int>(xprng.uniforme01() * kTercio12);
    if (xprng.uniforme01() < 0.5) std::swap(ind1, ind2);

    parent1 = &pool->parents[ind1];
    parent2 = &pool->parents[ind2];

    if (xprng.uniforme01() < pool->mutParams.frecRec) {
      const double r = xprng.uniforme01();
      int topePosInicio;
      if      (r < 0.1) topePosInicio = pool->poolParams.topePosInicio[0];
      else if (r < 0.5) topePosInicio = pool->poolParams.topePosInicio[1];
      else if (r < 0.9) topePosInicio = pool->poolParams.topePosInicio[2];
      else              topePosInicio = pool->poolParams.topePosInicio[3];

#if FLAG_POR_BLOQUES
      do {
        posiblock = static_cast<int>(xprng.uniforme01() * (parent1->nSeg - 1));
        posigen = parent1->segBl[posiblock] +
                  static_cast<int>(xprng.uniforme01() *
                                   (parent1->segBl[posiblock + 1] -
                                    parent1->segBl[posiblock]));
      } while (posigen > topePosInicio);
#else
      posigen =
          kResolucion + static_cast<int>(xprng.uniforme01() * topePosInicio);
#endif
      // 50%: intra-parent — use parent1 for both halves.
      if (xprng.uniforme01() < 0.5) parent2 = parent1;

      if (!CheckParentMutations(parent1, posigen)) continue;
      if (!CheckParentMutations(parent2, posigen)) continue;

      int nSeg = 0;
      // Inherit segments before posigen from parent1.
      for (int j = 0; j < parent1->nSeg; ++j) {
        if (parent1->segBl[j] < posigen) {
          child.segBl[nSeg] = parent1->segBl[j];
          child.NeBl[nSeg]  = parent1->NeBl[j];
          ++nSeg;
        } else {
          break;
        }
      }
      // Splice in the parent2 segment that contains posigen.
      for (int j = 1; j < parent2->nSeg + 1; ++j) {
        if (parent2->segBl[j] > posigen) {
          child.segBl[nSeg] = posigen;
          child.NeBl[nSeg]  = parent2->NeBl[j - 1];
          ++nSeg;
          break;
        }
      }
      // Inherit remaining segments after posigen from parent2.
      for (int j = 1; j < parent2->nSeg; ++j) {
        if (parent2->segBl[j] > posigen) {
          child.segBl[nSeg] = parent2->segBl[j];
          child.NeBl[nSeg]  = parent2->NeBl[j];
          ++nSeg;
        }
      }
      child.nSeg        = nSeg;
      child.segBl[nSeg] = parent2->segBl[parent2->nSeg];
    } else {
      child = *parent1;
    }

    isValidChild = MutateFusion(&child, pool);
  }

  MutateNeRnd(&child, pool);
  MutateLateral(&child, pool);

  child.efval = parent1->efval;

  // If the children pool isn't full, append. Otherwise replace the
  // current worst child if this one is better.
  if (childIdx < kNHijos) {
    child.SCval = CalculaSCScoreOnly(&child, sampleInfo);
    pool->children[childIdx] = child;
  } else {
    if (!worstChild->initialized) {
      RefreshWorstChild(pool, worstChild);
    }
    const double maxSC = worstChild->score;
    const int idxMaxSC = worstChild->index;
    child.SCval = CalculaSCScoreCutoff(&child, sampleInfo, maxSC);
    if (child.SCval < maxSC) {
      pool->children[idxMaxSC] = child;
      RefreshWorstChild(pool, worstChild);
    }
  }
}

// Lateral mutation: move each interior segment boundary by a small
// integer offset, biased small.
void MutateLateral(Bicho* bicho, Pool* pool) {
  for (int posiblock = 1; posiblock < bicho->nSeg - 1; ++posiblock) {
    int efecto;
    if (xprng.uniforme01() < pool->mutParams.frecMutLateral) {
      efecto =
          static_cast<int>(xprng.uniforme01() *
                           pool->mutParams.efectoMutLateral) +
          1;
    } else {
      efecto = 1;
      if (xprng.uniforme01() < 0.5) efecto = 0;
    }
    if (xprng.uniforme01() < 0.5) {
      if ((bicho->segBl[posiblock] - bicho->segBl[posiblock - 1]) >
          (kResolucion + efecto)) {
        bicho->segBl[posiblock] -= efecto;
      }
    } else {
      if ((bicho->segBl[posiblock + 1] - bicho->segBl[posiblock]) >
          (kResolucion + efecto)) {
        bicho->segBl[posiblock] += efecto;
      }
    }
  }
}

// Random Ne mutation: perturb each segment's Ne and reject if the
// resulting jump to either neighbour exceeds the allowed factor.
void MutateNeRnd(Bicho* bicho, Pool* pool) {
  for (int posiblock = 0; posiblock < bicho->nSeg; ++posiblock) {
    double efecto;
    if (xprng.uniforme01() < pool->mutParams.frecMut) {
      efecto = bicho->NeBl[posiblock] * xprng.uniforme01() *
               pool->mutParams.efectoMut;
    } else {
      efecto =
          bicho->NeBl[posiblock] * xprng.uniforme01() * kEfectoMutSuave;
    }
    if (xprng.uniforme01() < 0.5) efecto = -efecto;

    const double NeMutado = bicho->NeBl[posiblock] + efecto;

    if (posiblock > 0 &&
        NeMutado / bicho->NeBl[posiblock - 1] < kInvTopeSalto2) {
      continue;
    }
    if (posiblock < (bicho->nSeg - 1) &&
        NeMutado / bicho->NeBl[posiblock + 1] > kTopeSalto2) {
      continue;
    }

    bicho->NeBl[posiblock] =
        std::max(std::min(NeMutado, MAX_NE_SIZE), MIN_NE_SIZE);
  }
  // Clamp the last Ne to within ±20% of the previous one.
  const int lastIdx = bicho->nSeg - 1;
  const int prevToLastIdx = bicho->nSeg - 2;
  if (bicho->NeBl[lastIdx] > bicho->NeBl[prevToLastIdx] * 1.2) {
    bicho->NeBl[lastIdx] = bicho->NeBl[prevToLastIdx] * 1.2;
  } else if (bicho->NeBl[lastIdx] < bicho->NeBl[prevToLastIdx] / 1.2) {
    bicho->NeBl[lastIdx] = bicho->NeBl[prevToLastIdx] / 1.2;
  }
}

void ReproducePool(Pool* pool, GsampleInfo* sampleInfo) {
  // Produce the next generation of children.
  ChildWorstState worstChild = {-1, 0, false};
  for (int des = 0; des < pool->numChildren; ++des) {
    Reproduce(pool, sampleInfo, des, &worstChild);
  }

  // Promote children into the second third of the parent pool when
  // they beat one of those parents. The top third (kTercio1) of
  // parents is preserved untouched.
  if (pool->mutParams.solape) {
    for (int i = kTercio1; i < kTercio12; ++i) {
      double minSC = MAX_DOUBLE;
      int idxMinSC = 0;
      for (int j = 0; j < kNHijos; ++j) {
        if (pool->children[j].SCval < minSC) {
          minSC = pool->children[j].SCval;
          idxMinSC = j;
        }
      }
      bool foundBetterChild = false;
      for (int j = kTercio1; j < kTercio12; ++j) {
        if (minSC < pool->parents[j].SCval) {
          pool->parents[j] = pool->children[idxMinSC];
          pool->children[idxMinSC].SCval = MAX_DOUBLE;
          foundBetterChild = true;
          break;
        }
      }
      if (!foundBetterChild) break;
    }
  } else {
    // Wholesale replace parents with the first kTercio12 children.
    for (int i = 0; i < kTercio12; ++i) {
      pool->parents[i] = pool->children[i];
    }
  }

  // Inversion: swap the two trailing Ne values for any parent where
  // their ratio against the previous one stays within the jump factor.
  if (xprng.uniforme01() < kFrecInversion) {
    for (int i = 0; i < kTercio12; ++i) {
      Bicho* bicho = &pool->parents[i];
      if (bicho->nSeg >= 4) {
        const int posiblock = bicho->nSeg - 2;
        const double NeRate =
            bicho->NeBl[posiblock + 1] / bicho->NeBl[posiblock - 1];
        if (NeRate > kInvTopeSalto && NeRate < kTopeSalto &&
            bicho->NeBl[posiblock + 1] != bicho->NeBl[posiblock]) {
          std::swap(bicho->NeBl[posiblock + 1], bicho->NeBl[posiblock]);
          bicho->SCval = CalculaSCScoreOnly(bicho, sampleInfo);
        }
      }
    }
  }

  // Selection sort parents by ascending SCval. (This is a textbook
  // selection sort — not the quirky one in PrePopulatePool.)
  for (int i = 0; i < kTercio12 - 1; ++i) {
    int idxMinSC = i;
    for (int j = i + 1; j < kTercio12; ++j) {
      if (pool->parents[j].SCval < pool->parents[idxMinSC].SCval) {
        idxMinSC = j;
      }
    }
    if (idxMinSC != i) {
      std::swap(pool->parents[i], pool->parents[idxMinSC]);
    }
  }
}

// Mutation-rate schedule applied at specific GA generation milestones.
static void ApplyGenerationSchedule(Pool* pool, int gen) {
  switch (gen) {
    case 300:
      SetMutationParameters(pool, 0.2, 0.2, 2.0, 0.2, 0.5, 0.2,
                            pool->mutParams.maxSegmentos + 2, true);
      break;
    case 600:
      SetMutationParameters(pool, 0.5, 0.05, 1.0, 0.2, 0.5, 0.2,
                            pool->mutParams.maxSegmentos + 10, true);
      break;
    case 700:
      SetMutationParameters(pool, 0.5, 0.2, 1.0, 1.0, 0.0, 0.2,
                            pool->mutParams.maxSegmentos + 10, true);
      break;
    case 710:
      SetMutationParameters(pool, 0.5, 0.2, 1.0, 1.0, 0.0, 0.2,
                            pool->mutParams.maxSegmentos + 10, false);
      break;
    case 720:
      SetMutationParameters(pool, 0.5, 0.1, 1.0, 0.95, 0.0, 0.2,
                            pool->mutParams.maxSegmentos + 10, false);
      break;
    case 730:
      SetMutationParameters(pool, 0.5, 0.04, 1.0, 0.95, 0.0, 0.2,
                            pool->mutParams.maxSegmentos + 20, false);
      break;
    default:
      break;
  }
}

void RunDbg(Pool* pool, GsampleInfo* sampleInfo, std::string fichSal) {
  std::ofstream salida;
  const std::string fichero_sal_evol = fichSal + "_GONE2_evol";
  const std::string fichero_dbg      = fichSal + "_GONE2_dbg";

  salida.open(fichero_dbg, std::ios::app);
  for (int j = 0; j < pool->parents[0].nSeg; ++j) {
    salida << j << "\t" << pool->parents[0].segBl[j] << "\t"
           << pool->parents[0].NeBl[j] / 2.0 << "\n";
  }
  salida << "\n";
  salida.close();

  for (int gen = 0; gen < pool->numGenerations; ++gen) {
    // Original predicate `gen == int(gen/100.0)` is only true for
    // gen=0, so the dump only fires for the first generation.
    if (gen == static_cast<int>(static_cast<double>(gen) / 100.0)) {
      salida.open(fichero_dbg, std::ios::app);
      for (int j = 0; j < pool->parents[0].nSeg; ++j) {
        salida << "-> " << gen << "\t" << j << "\t"
               << pool->parents[0].segBl[j] << "\t"
               << pool->parents[0].NeBl[j] / 2.0 << "\n";
      }
      salida << "\n";
      salida.close();
    }
    ApplyGenerationSchedule(pool, gen);
    ReproducePool(pool, sampleInfo);

    double SCmed = 0, nSegMed = 0;
    for (int i = 0; i < kTercio12; ++i) {
      SCmed   += pool->parents[i].SCval;
      nSegMed += pool->parents[i].nSeg;
    }
    SCmed   /= kTercio12;
    nSegMed /= kTercio12;

    salida.open(fichero_sal_evol, std::ios::app);
    salida << gen << "\t" << pool->parents[0].SCval << "\t" << SCmed << "\t"
           << pool->parents[0].nSeg << "\t" << nSegMed << "\n";
    salida.close();
  }

  salida.open(fichero_dbg, std::ios::app);
  for (int j = 0; j < pool->parents[0].nSeg; ++j) {
    salida << j << "\t" << pool->parents[0].segBl[j] << "\t"
           << pool->parents[0].NeBl[j] / 2.0 << "\n";
  }
  salida.close();
}

void Run(Pool* pool, GsampleInfo* sampleInfo, std::string /*fichSal*/) {
  for (int gen = 0; gen < pool->numGenerations; ++gen) {
    ApplyGenerationSchedule(pool, gen);
    ReproducePool(pool, sampleInfo);
  }
}

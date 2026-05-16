// GONE2 — Genetic algorithm driver
// Wraps the per-round GA work into a parallel loop and writes the
// final averaged Ne / d² output.

#include "libgone.hpp"

#include <algorithm>
#include <limits>
#include <vector>

void gone(AppParams *params, std::string fichero, int argc, char *argv[],
          PopulationInfo *popInfo, SampleInfo *sInfo) {
  const double clow = 0, chigh = 0.5;

  if (params->hc <= params->lc) {
    std::cerr << " Invalid range of recombination frequencies." << std::endl;
    exit(1);
  }

  xprng.setSeed(params->semilla);

  const std::string fichsal = fichero.substr(0, fichero.length() - 7);
  const std::string fichero_sal_NeH  = fichsal + "_GONE2_Ne";
  const std::string fichero_sal_d2   = fichsal + "_GONE2_d2";
  const std::string fichero_sal_evol = fichsal + "_GONE2_evol";
  const std::string fichero_dbg      = fichsal + "_GONE2_dbg";
  std::ofstream salida;

  GsampleInfo *sampleInfo = new GsampleInfo();
  sampleInfo->flags     = params->flags;
  sampleInfo->hetEspAll = popInfo->hetEspAll;
  sampleInfo->hetEsp    = popInfo->hetEsp;
  sampleInfo->hetAvgAll = popInfo->hetAvgAll;
  sampleInfo->hetAvg    = popInfo->hetAvg;
  // Defaults read at compile time from -DGA_* flags; the combo path
  // (gone-ncurses_combo) overrides this between successive GA passes.
  sampleInfo->gaConfig  = MakeDefaultGAConfig();

  // -m: heterozygosity anchor. With per-site mutation rate μ and
  // genome-wide expected heterozygosity H = E[2pq] averaged over ALL
  // sites (segregating + invariant), mutation-drift balance gives
  // E[H] = θ/(1+θ) ≈ θ for small θ, where θ = 4·Ne_diploid·μ. We
  // anchor Necons (haploid units, = 2·Ne_diploid) toward
  // H/(2μ). The anchor is a soft squared-log penalty so an order-
  // of-magnitude mismatch in μ shifts the deep-end estimate but
  // doesn't lock it. Skipped for haploid data (haplotype == 1) where
  // 2pq isn't the right summary statistic.
  if (params->mutationRate > 0.0 && params->haplotype != 1 &&
      sampleInfo->hetEspAll > 0.0) {
    const double Ne_haploid_anchor =
        sampleInfo->hetEspAll / (2.0 * params->mutationRate);
    sampleInfo->gaConfig.anchorNeHaploid = Ne_haploid_anchor;
    sampleInfo->gaConfig.anchorLambda    = 1e-4;
    std::cerr << " Heterozygosity anchor: H=" << sampleInfo->hetEspAll
              << " μ=" << params->mutationRate
              << " → Ne_diploid≈" << Ne_haploid_anchor / 2.0
              << " (soft penalty λ=" << sampleInfo->gaConfig.anchorLambda
              << ")\n";
  }

  if (params->basecallerror == 0) {
    sampleInfo->basecallcorrec = 1;
  } else {
    sampleInfo->basecallcorrec =
        pow((1 - 4 * params->basecallerror * (1 - params->basecallerror)), 2);
  }
  popInfo->basecallcorrec = sampleInfo->basecallcorrec;

  sampleInfo->mix = sInfo->mix;
  if ((sampleInfo->flags & FLAG_RESIZE_BINS) > 0) {
    sampleInfo->sizeBins = params->sizeBins;
  } else {
    sampleInfo->sizeBins = DEFAULT_BIN_SIZE;
  }
  sampleInfo->nBins = kNumBins;

  const int linesRead = ProcessFile(fichero, clow, chigh, sampleInfo);
  if (linesRead < 15) {
    std::cerr
        << "There are not enough recombination bins to perform the analysis.\n";
    exit(1);
  }
  // Copy bin observations needed for the mix calculation.
  memcpy(&(sInfo->xc[0]), &(sampleInfo->cVal[0]), linesRead * sizeof(double));
  for (int i = 0; i < linesRead; ++i) {
    sInfo->nxc[i] = static_cast<long int>(sampleInfo->nBin[i]);
  }
  memcpy(&(sInfo->d2[0]), &(sampleInfo->d2cObs[0]), linesRead * sizeof(double));
  sInfo->binMax = linesRead;

  sampleInfo->cValMin = MAX_DOUBLE;

  if (sampleInfo->haplotype == 0) {
    sampleInfo->correccion = 1.0 / Square<double>(1.0 + sampleInfo->fVal);
  }

  sampleInfo->muestraSalida = params->muestraSalida;

  CalculateSumNBins(sampleInfo);
  CalculateAverageNe(sampleInfo);

  // Successive-generations sampling correction per bin.
  const double g = params->ngensampling;
  for (int i = 0; i < sampleInfo->nBins; ++i) {
    const int indx = sampleInfo->indx[i];
    const double ac = 1 - sampleInfo->cVal[indx];
    if (g == 1) {
      sampleInfo->ngensamplingcorrec[indx] = 1;
    } else {
      sampleInfo->ngensamplingcorrec[indx] =
          (g - 2 * ac + 2 * pow(ac, g + 1) - ac * ac * g) /
          (g * g * (ac - 1) * (ac - 1));
    }
  }

  if ((sampleInfo->flags & FLAG_DEBUG) > 0) {
    salida.open(fichero_sal_evol, std::ios::out);
    salida << "Gener\tSCbest\tSCmed1\tnsegbest\tnsegmed\n";
    salida.close();
  }

  // Single-population panmictic GA path. Mix is handled by the caller.
  if (!params->mix) {
    Pool *pool = new Pool();
    SetInitialPoolParameters(pool, sampleInfo->cValMin);
    pool->gaConfig = sampleInfo->gaConfig;
    PrePopulatePool(pool, sampleInfo);

    if ((sampleInfo->flags & FLAG_DEBUG) > 0) {
      salida.open(fichero_dbg, std::ios::app);
      for (int j = 0; j < pool->parents[0].nSeg; ++j) {
        salida << j << "\t" << pool->parents[0].segBl[j] << "\t"
               << pool->parents[0].NeBl[j] / 2.0 << "\n";
      }
      salida << "\n";
      salida.close();
    }

    // Configs to run. Non-combo builds run a single pass with the
    // compile-time-derived config already in sampleInfo->gaConfig
    // (so behaviour matches the pre-runtime-config code exactly).
    // The combo target tries trunc05_kick, L2, and L1+kicks in turn
    // and keeps the one whose averaged predicted d² fits best.
    std::vector<GAConfig> configs;
#ifdef GONE_COMBO
    configs.push_back(MakeComboTruncKickConfig());
    configs.push_back(MakeComboL2Config());
    configs.push_back(MakeComboL1KickConfig());
    // The heterozygosity anchor (if any) applies equally to all
    // three passes — we want the same deep-end target regardless of
    // which smoothness penalty is on.
    for (GAConfig& c : configs) {
      c.anchorNeHaploid = sampleInfo->gaConfig.anchorNeHaploid;
      c.anchorLambda    = sampleInfo->gaConfig.anchorLambda;
    }
#else
    configs.push_back(sampleInfo->gaConfig);
#endif

    double avgD2Pred[kNumLinMax] = {};
    double avgNe[MAXBINS] = {};
    // Best-so-far across passes (only one pass in non-combo builds).
    double bestAvgD2Pred[kNumLinMax] = {};
    double bestAvgNe[MAXBINS] = {};
    int    bestGmax2 = 0;
    double bestResidual = std::numeric_limits<double>::infinity();
    std::string bestLabel;

    static const Pool emptyPool = Pool();
    int numThreads = params->numThreads;
    if (numThreads == 0) {
      numThreads = omp_get_max_threads();
    }

    // Per-thread accumulators.
    //   tavgNe[t][i]      — indexed per generation, up to conta2 which can
    //                       reach gmax ≤ kNumGenMax-10. Reduced by
    //                       geometric mean, so start filled with 1.0.
    //   tavgD2Pred[t][i]  — indexed per recombination bin, i < nBins ≤ 60.
    double **tavgNe     = new double *[numThreads];
    double **tavgD2Pred = new double *[numThreads];
    int    *maxNeConta  = new int[numThreads]{};
    for (int z = 0; z < numThreads; ++z) {
      tavgD2Pred[z] = new double[kNumLinMax]{};
      tavgNe[z]     = new double[kNumGenMax]{};
      std::fill_n(tavgNe[z], kNumGenMax, 1.0);
    }
    params->progress.InitCurrentTask(
        static_cast<float>(GONE_ROUNDS * configs.size()));
    std::vector<double> liveLogNe(MAXBINS, 0.0);
    std::vector<int> liveNeCounts(MAXBINS, 0);
    int completedRounds = 0;
    int liveMaxConta = 0;
    const std::string neUnitLabel =
        params->haplotype == 1 ? "Ne_haploids" : "Ne_diploids";

    for (size_t passIdx = 0; passIdx < configs.size(); ++passIdx) {
      const GAConfig& cfg = configs[passIdx];
      sampleInfo->gaConfig = cfg;
      // Reseed the global RNG at the start of every pass so each pass
      // sees the same random sequence as its standalone binary would.
      // Without this, pass 2 / pass 3 inherit whatever state pass 1
      // left, the GAs land in different basins than the individual
      // builds, and the curves pick up the "sudden jumps" the user
      // reported.
      xprng.setSeed(params->semilla);
      // Reset per-pass accumulators so the GA starts fresh under
      // this config.
      std::fill_n(avgD2Pred, kNumLinMax, 0.0);
      std::fill_n(avgNe, MAXBINS, 0.0);
      for (int z = 0; z < numThreads; ++z) {
        std::fill_n(tavgD2Pred[z], kNumLinMax, 0.0);
        std::fill_n(tavgNe[z], kNumGenMax, 1.0);
        maxNeConta[z] = 0;
      }
      std::fill(liveLogNe.begin(), liveLogNe.end(), 0.0);
      std::fill(liveNeCounts.begin(), liveNeCounts.end(), 0);
      int passCompleted = 0;
      liveMaxConta = 0;
      // The score box still tracks the per-round running min of the
      // un-penalized bin-weighted d² residual (no smoothness term, so
      // values stay comparable across algorithms during the live run).
      // The combo picker, by contrast, uses the residual of each
      // pass's averaged Ne curve — computed once after the round
      // loop — so the chosen pass matches what gets written to
      // _GONE2_Ne. See the block below the parallel section.
      double passBestResidual = std::numeric_limits<double>::infinity();
      params->progress.ResetBestScore();
      const std::string passPrefix =
          configs.size() > 1
              ? "[" + std::to_string(passIdx + 1) + "/" +
                std::to_string(configs.size()) + " " +
                std::string(cfg.label) + "] "
              : std::string();

#pragma omp parallel
    {
      const int tid = omp_get_thread_num();
      Pool *privpool = new Pool();

      for (int _i = tid; _i < GONE_ROUNDS; _i += numThreads) {
        *privpool = emptyPool;
        SetInitialPoolParameters(privpool, sampleInfo->cValMin);
        // SetInitialPoolParameters resets pool->gaConfig to the
        // compile-time defaults; in the combo build that wipes the
        // per-pass kick/scout values we just put on sampleInfo. Copy
        // them back so MutateNeRnd / Run actually see this pass's
        // configuration rather than the baseline.
        privpool->gaConfig = sampleInfo->gaConfig;
        PrePopulatePool(privpool, sampleInfo);
        if ((sampleInfo->flags & FLAG_DEBUG) > 0) {
          RunDbg(privpool, sampleInfo, fichsal);
        } else {
          Run(privpool, sampleInfo, fichsal);
        }

        Bicho *bestBicho = &privpool->parents[0];
        double bestD2Pred[kMaxD2PredBins] = {};
        CalculaSC(bestBicho, sampleInfo, bestD2Pred);
        for (int i = 0; i < sampleInfo->nBins; ++i) {
          tavgD2Pred[tid][i] += bestD2Pred[i] / GONE_ROUNDS;
        }

#ifdef GA_LOG_ROUNDS
        // Append one line to fichsal + "_GA_rounds.csv":
        //   round,scval,nSeg,ne_at_gen1,ne_at_gen10,ne_at_gen100
        {
          const std::string fichRounds = fichsal + "_GA_rounds.csv";
#pragma omp critical(ga_log_rounds)
          {
            static bool headerWritten = false;
            std::ofstream out(fichRounds, std::ios::app);
            if (!headerWritten) {
              out << "round,scval,nSeg,ne_gen1,ne_gen10,ne_gen100\n";
              headerWritten = true;
            }
            // Look up Ne at three sentinel generations from the best
            // bicho's piecewise-constant schedule.
            auto neAtGen = [&](int gen) -> double {
              int seg = 0;
              while (seg < bestBicho->nSeg &&
                     bestBicho->segBl[seg + 1] <= gen) ++seg;
              if (seg >= bestBicho->nSeg) seg = bestBicho->nSeg - 1;
              const double scale = params->haplotype == 1 ? 1.0 : 0.5;
              return std::max(MIN_NE_SIZE, bestBicho->NeBl[seg] * scale);
            };
            out << (_i + 1) << "," << bestBicho->SCval << ","
                << bestBicho->nSeg << "," << neAtGen(1) << ","
                << neAtGen(10) << "," << neAtGen(100) << "\n";
          }
        }
#endif

        // Count of generations spanned by the best bicho's segments.
        int conta = 0;
        for (int i = 0; i < bestBicho->nSeg; ++i) {
          for (int j = bestBicho->segBl[i]; j < bestBicho->segBl[i + 1]; ++j) {
            ++conta;
          }
        }
        if (conta > kNumLinMax) conta = kNumLinMax;
        const int conta2 = conta;
        double sumNe[kNumGenMax] = {0};

        // Geometric mean over the top muestraSalida parents.
        for (int i = 0; i < conta; ++i) sumNe[i] = 1;
        for (int ii = 0; ii < sampleInfo->muestraSalida; ++ii) {
          conta = 0;
          for (int i = 0; i < privpool->parents[ii].nSeg; ++i) {
            for (int j = privpool->parents[ii].segBl[i];
                 j < privpool->parents[ii].segBl[i + 1]; ++j) {
              sumNe[conta] *= pow(privpool->parents[ii].NeBl[i],
                                  1.0 / sampleInfo->muestraSalida);
              ++conta;
              if (conta > privpool->poolParams.gmax[0]) break;
            }
            if (conta > privpool->poolParams.gmax[0]) break;
          }
        }

        // Convert to diploid units and accumulate geometric mean across rounds.
        for (int i = 0; i < conta2; ++i) {
          sumNe[i] /= 2;
          tavgNe[tid][i] *= pow(sumNe[i], 1.0 / GONE_ROUNDS);
        }
        std::vector<double> roundNe(conta2, MIN_NE_SIZE);
        for (int i = 0; i < conta2; ++i) {
          const double reportedNe =
              params->haplotype == 1 ? sumNe[i] * 2.0 : sumNe[i];
          roundNe[i] = std::max(reportedNe, MIN_NE_SIZE);
        }
        if (conta2 > maxNeConta[tid]) maxNeConta[tid] = conta2;
#pragma omp critical(gone_progress_update)
        {
          ++completedRounds;
          ++passCompleted;
          liveMaxConta = std::max(liveMaxConta, conta2);
          for (int i = 0; i < conta2; ++i) {
            liveLogNe[i] += std::log(roundNe[i]);
            liveNeCounts[i] += 1;
          }
          std::vector<double> liveNe(liveMaxConta, MIN_NE_SIZE);
          for (int i = 0; i < liveMaxConta; ++i) {
            if (liveNeCounts[i] > 0) {
              liveNe[i] = std::exp(liveLogNe[i] / liveNeCounts[i]);
            }
          }
          // Bin-weighted d² residual for THIS round's best bicho,
          // with the smoothness penalty stripped (the score box
          // would otherwise be incomparable across passes that use
          // different penalties).
          double roundRes = 0.0, roundN = 0.0;
          for (int i = 0; i < sampleInfo->nBins; ++i) {
            if (sampleInfo->cVal[i] != 0) {
              const double d = sampleInfo->d2cObs[i] - bestD2Pred[i];
              roundRes += sampleInfo->nBin[i] * d * d;
              roundN   += sampleInfo->nBin[i];
            }
          }
          const double roundScore =
              roundN > 0
                  ? roundRes / roundN
                  : std::numeric_limits<double>::infinity();
          if (roundScore < passBestResidual) {
            passBestResidual = roundScore;
          }
          params->progress.SetStatusDetail(
              passPrefix + "Completed GA round " +
              std::to_string(passCompleted) + " of " +
              std::to_string(GONE_ROUNDS));
          params->progress.SetNeSnapshot(liveNe, neUnitLabel);
          params->progress.SetBestScore(roundScore, passCompleted,
                                        GONE_ROUNDS);
          params->progress.SetTaskProgress(
              static_cast<float>(completedRounds));
        }
      }
      delete privpool;
    }

    // Reduction across threads for this pass.
    std::fill_n(avgNe, MAXBINS, 1.0);
    for (int nt = 0; nt < numThreads; ++nt) {
      for (int i = 0; i < maxNeConta[nt]; ++i) {
        avgNe[i] *= tavgNe[nt][i];
      }
      for (int i = 0; i < sampleInfo->nBins; ++i) {
        avgD2Pred[i] += tavgD2Pred[nt][i];
      }
    }

    // Score this pass by the residual of the averaged Ne curve's
    // own d² prediction — the curve that will be saved to
    // _GONE2_Ne. Build a synthetic Bicho with one segment per
    // generation from avgNe[] (diploid units → ×2 to haploid for
    // Bicho::NeBl), then call CalculaSC to get its predicted d².
    // This replaces the older running-min-of-per-round-residuals
    // criterion, which could pick a pass whose averaged curve fit
    // worse than another's (the per-round best was a lucky outlier).
    {
      int curveLen = 0;
      for (int nt = 0; nt < numThreads; ++nt) {
        if (maxNeConta[nt] > curveLen) curveLen = maxNeConta[nt];
      }
      const int maxSeg = std::min(curveLen, kNumLinMax - 1);
      if (maxSeg >= 1) {
        Bicho synth = {};
        synth.nSeg = maxSeg;
        for (int i = 0; i < maxSeg; ++i) {
          synth.segBl[i] = i;
          synth.NeBl[i] = std::max(MIN_NE_SIZE, avgNe[i] * 2.0);
        }
        synth.segBl[maxSeg] = maxSeg;
        synth.efval = sampleInfo->fVal;
        double predD2[kMaxD2PredBins] = {};
        CalculaSC(&synth, sampleInfo, predD2);
        double res = 0.0, n = 0.0;
        for (int i = 0; i < sampleInfo->nBins; ++i) {
          if (sampleInfo->cVal[i] != 0) {
            const double d = sampleInfo->d2cObs[i] - predD2[i];
            res += sampleInfo->nBin[i] * d * d;
            n   += sampleInfo->nBin[i];
          }
        }
        passBestResidual =
            (n > 0) ? res / n : std::numeric_limits<double>::infinity();
      }
    }

    if (passBestResidual < bestResidual) {
      bestResidual = passBestResidual;
      bestLabel    = cfg.label ? cfg.label : "";
      bestGmax2    = pool->poolParams.gmax[2];
      std::copy(avgNe, avgNe + MAXBINS, bestAvgNe);
      std::copy(avgD2Pred, avgD2Pred + kNumLinMax, bestAvgD2Pred);
    }
    }  // end for passIdx

    // Free per-thread accumulators after all passes are done.
    for (int nt = 0; nt < numThreads; ++nt) {
      delete[] tavgNe[nt];
      delete[] tavgD2Pred[nt];
    }
    delete[] tavgNe;
    delete[] tavgD2Pred;
    delete[] maxNeConta;

    // Push the winning curve to the TUI so the user sees the chosen
    // result after the program finishes. In non-combo builds this is
    // a no-op (the last in-loop snapshot already carries the same
    // values); in combo it might roll back from a later, worse pass.
    // The score box also rolls back to the winning pass's residual
    // so the displayed number matches the kept solution rather than
    // sitting at whatever the last-running pass produced.
    if (configs.size() > 1) {
      std::vector<double> finalNe(bestGmax2 + 1, MIN_NE_SIZE);
      for (int i = 0; i <= bestGmax2 && i < kNumGenMax; ++i) {
        const double scale = params->haplotype == 1 ? 2.0 : 1.0;
        finalNe[i] = std::max(MIN_NE_SIZE, bestAvgNe[i] * scale);
      }
      params->progress.SetNeSnapshot(finalNe, neUnitLabel);
      params->progress.SetStatusDetail(
          "Best fit: " + bestLabel +
          " (bin-weighted d² residual " + std::to_string(bestResidual) + ")");
      params->progress.SetBestScoreLabeled(
          bestResidual, "Best fit", "kept: " + bestLabel);
      std::cerr << " Combo kept: " << bestLabel
                << " (bin-weighted d² residual "
                << std::to_string(bestResidual) << ")\n";
    }

    // Write the Ne output of the WINNING pass. Internal Ne is haploid;
    // we divide by 2 to report diploid Ne (except for the haploid case
    // which multiplies back by 2 here).
    salida.open(fichero_sal_NeH, std::ios::out);
    if (params->haplotype == 1) {
      salida << "Generation\tNe_haploids\n";
      for (int j = 0; j < bestGmax2 + 1; ++j) {
        const int generacion = j + 1;
        if (generacion < 151) {
          bestAvgNe[j] *= 2;
          salida << generacion << "\t" << std::max(bestAvgNe[j], MIN_NE_SIZE)
                 << "\n";
        }
      }
    } else {
      salida << "Generation\tNe_diploids\n";
      for (int j = 0; j < bestGmax2 + 1; ++j) {
        const int generacion = j + 1;
        if (generacion < 151) {
          salida << generacion << "\t" << std::max(bestAvgNe[j], MIN_NE_SIZE)
                 << "\n";
        }
      }
    }
    salida.close();

    // Write d² (observed vs winning-pass predicted) per bin.
    salida.open(fichero_sal_d2, std::ios::out);
    salida << "c_bin\tnumber_of_SNP_pairs\tObserved_d2\tPredicted_d2\n";
    for (int i = 0; i < linesRead; ++i) {
      if (sampleInfo->cVal[i] != 0) {
        salida << std::fixed << std::setprecision(8) << sampleInfo->cVal[i]
               << "\t" << std::fixed << std::setprecision(0)
               << sampleInfo->nBin[i] << "\t" << std::fixed
               << std::setprecision(8) << sampleInfo->d2cObs[i] << "\t"
               << bestAvgD2Pred[i] << "\n";
      }
    }
    salida.close();

    delete pool;
    delete sampleInfo;
  }
}

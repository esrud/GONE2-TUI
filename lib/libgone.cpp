// GONE2 — Genetic algorithm driver
// Wraps the per-round GA work into a parallel loop and writes the
// final averaged Ne / d² output.

#include "libgone.hpp"

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

    double avgD2Pred[kNumLinMax] = {};
    double avgNe[MAXBINS] = {};

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
    params->progress.InitCurrentTask(static_cast<float>(GONE_ROUNDS));
    std::vector<double> liveLogNe(MAXBINS, 0.0);
    std::vector<int> liveNeCounts(MAXBINS, 0);
    int completedRounds = 0;
    int liveMaxConta = 0;
    const std::string neUnitLabel =
        params->haplotype == 1 ? "Ne_haploids" : "Ne_diploids";

#pragma omp parallel
    {
      const int tid = omp_get_thread_num();
      Pool *privpool = new Pool();

      for (int _i = tid; _i < GONE_ROUNDS; _i += numThreads) {
        *privpool = emptyPool;
        SetInitialPoolParameters(privpool, sampleInfo->cValMin);
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
          params->progress.SetStatusDetail(
              "Completed GA round " + std::to_string(completedRounds) +
              " of " + std::to_string(GONE_ROUNDS));
          params->progress.SetNeSnapshot(liveNe, neUnitLabel);
          params->progress.SetBestScore(bestBicho->SCval, completedRounds,
                                        GONE_ROUNDS);
          params->progress.SetTaskProgress(
              static_cast<float>(completedRounds));
        }
      }
      delete privpool;
    }

    // Reduction across threads.
    std::fill_n(avgNe, MAXBINS, 1.0);
    for (int nt = 0; nt < numThreads; ++nt) {
      for (int i = 0; i < maxNeConta[nt]; ++i) {
        avgNe[i] *= tavgNe[nt][i];
      }
      for (int i = 0; i < sampleInfo->nBins; ++i) {
        avgD2Pred[i] += tavgD2Pred[nt][i];
      }
      delete[] tavgNe[nt];
      delete[] tavgD2Pred[nt];
    }
    delete[] tavgNe;
    delete[] tavgD2Pred;
    delete[] maxNeConta;

    // Write the Ne output. Internal Ne is haploid; we divide by 2 to
    // report diploid Ne (except for the haploid case which multiplies
    // back by 2 here).
    salida.open(fichero_sal_NeH, std::ios::out);
    if (params->haplotype == 1) {
      salida << "Generation\tNe_haploids\n";
      for (int j = 0; j < pool->poolParams.gmax[2] + 1; ++j) {
        const int generacion = j + 1;
        if (generacion < 151) {
          avgNe[j] *= 2;
          salida << generacion << "\t" << std::max(avgNe[j], MIN_NE_SIZE)
                 << "\n";
        }
      }
    } else {
      salida << "Generation\tNe_diploids\n";
      for (int j = 0; j < pool->poolParams.gmax[2] + 1; ++j) {
        const int generacion = j + 1;
        if (generacion < 151) {
          salida << generacion << "\t" << std::max(avgNe[j], MIN_NE_SIZE)
                 << "\n";
        }
      }
    }
    salida.close();

    // Write d² (observed vs predicted) per bin.
    salida.open(fichero_sal_d2, std::ios::out);
    salida << "c_bin\tnumber_of_SNP_pairs\tObserved_d2\tPredicted_d2\n";
    for (int i = 0; i < linesRead; ++i) {
      if (sampleInfo->cVal[i] != 0) {
        salida << std::fixed << std::setprecision(8) << sampleInfo->cVal[i]
               << "\t" << std::fixed << std::setprecision(0)
               << sampleInfo->nBin[i] << "\t" << std::fixed
               << std::setprecision(8) << sampleInfo->d2cObs[i] << "\t"
               << avgD2Pred[i] << "\n";
      }
    }
    salida.close();

    delete pool;
    delete sampleInfo;
  }
}


#include "sample.hpp"
#include "pool.hpp"
#include <cmath>
#include <iomanip>
#include <ios>
#include <string>
#include <string.h>
#include <vector>

static inline int PopCount64(uint64_t value) {
  return __builtin_popcountll(value);
}

// Compute the non-reference allele frequency and non-reference
// homozygote frequency at each segregating locus in the sample. Marks
// loci as segregating only when the allele frequency is strictly in
// (0, 1) and at most 20% of genotype calls at that locus are missing.
void CalculateFrequencies(PopulationInfo* popInfo, SampleInfo* sampleInfo) {
  double contaIndX = 0;
  int locus;
  int indi;
  double counter = 0;
  uint8_t ff, gg;

  sampleInfo->hetEsp = 0;
  // 20% missing-data cutoff: require at least 80% of 2× sampleSizeIndi
  // observed allele calls.
  const int mincontaIndX = int(float(sampleInfo->sampleSizeIndi) * 2.0 * 0.8);
  for (int j = 0; j < sampleInfo->numLoci; ++j) {
    locus = sampleInfo->shuffledLoci[j];
    sampleInfo->frec[locus]    = 0;
    sampleInfo->homo[locus]    = 0;
    sampleInfo->segrega[locus] = false;
    contaIndX = 0;
    if (popInfo->haplotype != 2) {
      // Diploid unphased / haploid / pseudohaploid: one row per indi.
      for (int i = 0; i < sampleInfo->sampleSizeIndi; i++) {
        indi = sampleInfo->shuffledIndi[i];
        ff = popInfo->indi[indi][locus];
        if (ff < 9) {
          sampleInfo->frec[locus] += static_cast<double>(ff);
          if (ff == 2) ++sampleInfo->homo[locus];
          ++(++contaIndX);
        }
      }
    } else {
      // Phased diploids: two adjacent rows per individual.
      for (int i = 0; i < sampleInfo->sampleSizeIndi; ++(++i)) {
        indi = sampleInfo->shuffledIndi[i];
        ff = popInfo->indi[indi][locus];
        indi = sampleInfo->shuffledIndi[i + 1];
        gg = popInfo->indi[indi][locus];
        if ((ff < 9) and (gg < 9)) {
          sampleInfo->frec[locus] += static_cast<double>(ff);
          ++(++contaIndX);
          sampleInfo->frec[locus] += static_cast<double>(gg);
          ++(++contaIndX);
          if ((ff == 2) && (gg == 2)) ++(++sampleInfo->homo[locus]);
        }
      }
    }

    if (contaIndX == 0) {
      std::cerr << "There is no genotyping data for at least one SNP\n";
      exit(EXIT_FAILURE);
    }
    // Segregating filter: non-fixed allele frequency and enough non-
    // missing calls to clear the 20% missingness cutoff.
    if ((sampleInfo->frec[locus] > 0) && (contaIndX >= mincontaIndX) &&
        (sampleInfo->frec[locus] < contaIndX)) {
      ++sampleInfo->numSegLoci;
      sampleInfo->segrega[locus] = true;
      sampleInfo->frec[locus] /= contaIndX;
      sampleInfo->homo[locus] /= (contaIndX / 2.0);
      sampleInfo->hetEsp +=
          2 * sampleInfo->frec[locus] * (1 - sampleInfo->frec[locus]);
    }
    popInfo->avgNumIndiAnalyzed += static_cast<double>(contaIndX);
    counter++;

    if (sampleInfo->numSegLoci >= sampleInfo->numLoci) break;
  }
  popInfo->avgNumIndiAnalyzed /= counter * 2.0;
  sampleInfo->hetEsp /= static_cast<double>(sampleInfo->numSegLoci);
  sampleInfo->hetEspAll = sampleInfo->hetEsp *
                          static_cast<double>(sampleInfo->numSegLoci) /
                          static_cast<double>(sampleInfo->numLoci);
  popInfo->hetEsp    = sampleInfo->hetEsp;
  popInfo->hetEspAll = sampleInfo->hetEspAll;
}

// Sample-level inbreeding coefficient f and the expected
// heterozygosity hetpq = mean(2pq) over segregating loci.
void CalculateF(SampleInfo* sampleInfo, PopulationInfo* popInfo) {
  (void)popInfo;
  double acupq = 0;
  double acuPp2 = 0;
  for (int j = 0; j < sampleInfo->numLoci; ++j) {
    const int locus = sampleInfo->shuffledLoci[j];
    if (sampleInfo->segrega[locus]) {
      acuPp2 +=
          sampleInfo->homo[locus] - Square<double>(sampleInfo->frec[locus]);
      acupq += sampleInfo->frec[locus] * (1.0 - sampleInfo->frec[locus]);
    }
  }
  sampleInfo->f     = acuPp2 / acupq;
  sampleInfo->hetpq = acupq / sampleInfo->numSegLoci;
}


bool CalculateHeterozigosity(PopulationInfo* popInfo, SampleInfo* sampleInfo) {
    // Distribucion de heterocigosidades de los individuos (ver
    // Medidas_de_Parent_y_Het.docx)
    int* pIndi = sampleInfo->shuffledIndi;
    int* pLocus;

    int contaLocX;
    int i, j;
    bool isSegLocus;
    uint8_t ff,gg;

    if (popInfo->haplotype ==0){
      sampleInfo->hetAvg = 0;
      for (i = 0; i < sampleInfo->sampleSizeIndi; ++i) {
        double hetIndi=0;
        pLocus = sampleInfo->shuffledLoci;
        contaLocX = 0;
        for (j = 0; j < sampleInfo->numSegLoci; ++j) {
          isSegLocus = false;
          while (!isSegLocus) {
            if (sampleInfo->segrega[*pLocus]) {
              ff = popInfo->indi[*pIndi][*pLocus];
              if (ff == 1) {
                // acumulador de heterocigosidad de cada indiv
                ++hetIndi;
                ++contaLocX;
              } else if (ff < 9) {
                ++contaLocX;
              }
              isSegLocus = true;
            }
            ++pLocus;
          }
        }
        if (contaLocX == 0) {
          std::cerr << "There is no genotyping data for at least one individual"
                    << std::endl;
          delete popInfo;
          delete sampleInfo;
          return false;
        }
        hetIndi /= contaLocX;  // het de cada individuo
        sampleInfo->hetAvg += hetIndi;
        ++pIndi;
      }
      sampleInfo->hetAvg /= static_cast<double>(sampleInfo->sampleSizeIndi);
      sampleInfo->hetAvgAll = sampleInfo->hetAvg * static_cast<double>(sampleInfo->numSegLoci) / static_cast<double>(sampleInfo->numLoci);
      popInfo->hetAvg=sampleInfo->hetAvg;
      popInfo->hetAvgAll = sampleInfo->hetAvgAll;
    }
    else{ // PHASED
      sampleInfo->hetAvg = 0;
      for (i = 0; i < sampleInfo->sampleSizeIndi; ++(++i)) {
        double hetIndi=0;
        pLocus = sampleInfo->shuffledLoci;
        contaLocX = 0;
        for (j = 0; j < sampleInfo->numSegLoci; ++j) {
          isSegLocus = false;
          while (!isSegLocus) {
            if (sampleInfo->segrega[*pLocus]) {
              ff = popInfo->indi[*pIndi][*pLocus];
              gg = popInfo->indi[*pIndi+1][*pLocus];
              if (ff+gg == 2) {
                // acumulador de heterocigosidad de cada indiv
                ++hetIndi;
                ++contaLocX;
              } else if ((ff < 9) && (gg < 9)) {
                ++contaLocX;
              }
              isSegLocus = true;
            }
            ++pLocus;
          }
        }
        if (contaLocX == 0) {
          std::cerr << "There is no genotyping data for at least one individual"
                    << std::endl;
          delete popInfo;
          delete sampleInfo;
          return false;
        }
        hetIndi /= contaLocX;  // het de cada individuo
        sampleInfo->hetAvg += hetIndi;
        ++(++pIndi);
      }
      sampleInfo->hetAvg /= (static_cast<double>(sampleInfo->sampleSizeIndi) / 2.0);
      sampleInfo->hetAvgAll = sampleInfo->hetAvg * static_cast<double>(sampleInfo->numSegLoci) / static_cast<double>(sampleInfo->numLoci);
      popInfo->hetAvg=sampleInfo->hetAvg;
      popInfo->hetAvgAll = sampleInfo->hetAvgAll;
    }

    return true;
}

namespace {

// Generations beyond gen 150 are not charted (matches the chart's
// kMaxGen in ncurses_tui.cpp and the published Ne file's gen cap).
constexpr int kChartMaxGen = 150;

// Build a naive single-population Ne curve (Hayes-style) from current
// per-bin d² accumulators and push it through ProgressStatus so the
// live chart has something to show while the slow d² loop runs.
//
// Called from inside the loop's progress-update critical section.
// Other threads may still be writing to the per-thread accumulators
// concurrently — the reduction is therefore a brief snapshot that
// can be slightly inconsistent across bins, but is good enough for
// visualisation. The published curve is overwritten by the GA's
// refined estimate (non-mix) or the Fst-corrected curve (mix) at
// the end of the run, so any transient inaccuracy here washes out.
void PushNaiveNeFromAccumulators(int numThreads,
                                 const int* binMaxes,
                                 long int** tnxc,
                                 double** txc,
                                 double** txD2,
                                 double** txW,
                                 ProgressStatus& progress,
                                 const std::string& unitLabel) {
  long int nxc[MAXBINS] = {0};
  double xc[MAXBINS] = {0}, xD2[MAXBINS] = {0}, xW[MAXBINS] = {0};
  int binMax = 0;
  for (int t = 0; t < numThreads; ++t) {
    if (binMaxes[t] > binMax) binMax = binMaxes[t];
    for (int b = 0; b <= binMaxes[t]; ++b) {
      nxc[b] += tnxc[t][b];
      xc[b]  += txc[t][b];
      xD2[b] += txD2[t][b];
      xW[b]  += txW[t][b];
    }
  }
  long int totalPairs = 0;
  for (int b = 0; b <= binMax; ++b) totalPairs += nxc[b];
  if (totalPairs < 50) return;
  std::vector<double> chartNe(kChartMaxGen, 0.0);
  for (int b = 0; b <= binMax; ++b) {
    if (nxc[b] == 0 || xW[b] <= 0) continue;
    const double cMean = 1.0 / (xc[b] / nxc[b]);  // harmonic mean of c
    const double d2 = xD2[b] / xW[b];
    if (d2 <= 0 || cMean <= 0 || cMean >= 0.5) continue;
    const double c2 = cMean * cMean;
    const double c12 = (1.0 - cMean) * (1.0 - cMean);
    const double NeRaw = (1.0 + c2 - 2.2 * d2 * c12) /
                         (2.0 * d2 * (1.0 - c12));
    if (NeRaw <= 0) continue;
    const int gen = static_cast<int>(std::round(1.0 / (2.0 * cMean))) - 1;
    if (gen >= 0 && gen < kChartMaxGen) {
      chartNe[gen] = std::max(1.0, NeRaw / 2.0);  // /2 for diploid report
    }
  }
  // Forward-fill so the curve renders continuously between sampled gens.
  double last = 0;
  for (int g = 0; g < kChartMaxGen; ++g) {
    if (chartNe[g] > 0) last = chartNe[g];
    else if (last > 0) chartNe[g] = last;
  }
  if (last > 0) progress.SetNeSnapshot(chartNe, unitLabel);
}

}  // namespace


void CalculateD2Parallel(PopulationInfo* popInfo, SampleInfo* sampleInfo,
                 AppParams* params) {
  // Calculo de D2:

  int* valid_idx = new int[popInfo->numLoci];
  int counter = 0;
  int i, j, k;
  int _i, _irepe;
  double gBase;
  double xlc = 0, xhc = 0;
  long int refextra;
  double infx, supx;
  int bin, extrabin;
  double epsilon=0;
  double iota=0;
  double kappa=0;
  double desde;
  int genot1,genot2;
  bool okchromsize;

  params->progress.InitCurrentTask(sampleInfo->numSegLoci - 1);

  // Inicializamos la tabla de índices válidos (en los que segrega[x] es true)
  for (int idx = 0; idx < popInfo->numLoci; idx++) {
    if (sampleInfo->segrega[sampleInfo->shuffledLoci[idx]] &&
        (sampleInfo->frec[sampleInfo->shuffledLoci[idx]] > params->MAF) &&
        ((1.0 - sampleInfo->frec[sampleInfo->shuffledLoci[idx]]) >
         params->MAF)) {
      valid_idx[counter] = sampleInfo->shuffledLoci[idx];
      counter++;
    }
  }

  gBase = 1.0 / (2.0 * params->hc);
  // gBase = std::max((1.0 / (2.0 * params->hc)), 1.0);
  //  binBase = trunc((gBase - 10.0) / 5.0) + 6.0;

  if (params->cMMb > 0) {
    for (i = 0; i < popInfo->numLoci; ++i) {
      popInfo->posiCM[i] = params->cMMb * popInfo->posiBP[i] / 1000000;
    }
    popInfo->Mtot = params->cMMb * popInfo->Mbtot / 100;
  }

  // Checking chomosome sizes 
  j=0;
  k=0;
  okchromsize=true;
  desde=popInfo->posiCM[j];
  for (i = 1; i < popInfo->numLoci; ++i) {
    if ((popInfo->cromo[i]!=popInfo->cromo[j]) || (i==(popInfo->numLoci-1))){
      if((popInfo->posiCM[i-1]-desde) < 20.0){
        std::cerr << "Chromosome "<< popInfo->cromo[j]<<" is too small. All chromosomes must be larger than 20cM" << std::endl;
        okchromsize=false;        
      }
      popInfo->chrsize[k] = popInfo->posiCM[i-1]-desde; // Se miden los cromosomas
      ++k;
      j=i;
      desde=popInfo->posiCM[j];
    }
  }
  if (!okchromsize){
    exit(EXIT_FAILURE);
  }

  if (params->distance == 0) {  // SIN TRANSFORMAR
    xlc = params->lc;
    xhc = params->hc;
  }
  if (params->distance == 1) {  // Haldane
    if (params->lc >= 0.5) {
      xlc = 100.0;
    } else {
      xlc = -log(1.0 - 2.0 * params->lc) / 2.0;
    }
    if (params->hc >= 0.5) {
      xhc = 200.0;
    } else {
      xhc = -log(1.0 - 2.0 * params->hc) / 2.0;
    }
  }
  if (params->distance == 2) {  // Kosambi
    if (params->lc >= 0.5) {
      xlc = 100.0;
    } else {
      xlc = -log((1.0 - 2.0 * params->lc) / (1.0 + 2.0 * params->lc)) / 4.0;
    }
    if (params->hc >= 0.5) {
      xhc = 200.0;
    } else {
      xhc = -log((1.0 - 2.0 * params->hc) / (1.0 + 2.0 * params->hc)) / 4.0;
    }
  }
  infx = xlc * 100.0;  // en cM
  supx = xhc * 100.0;  // en cM

  int locus1;
  int locus2;
  int* pIndi;

  double cenmor = 0;
  double contaIndX = 0;
  double tacuHoHo, tacuHoHetHetHo, tacuHetHet;
  uint8_t ss;

  double D, W,D2,W2;
  double recrate, generacion;

  int numThreads = params->numThreads;
  if (numThreads == 0) {
    numThreads = omp_get_max_threads();
  }

  double xW[MAXBINS] = {0}, xD2[MAXBINS] = {0};
  long int nxc[MAXBINS] = {0};
  double xc[MAXBINS] = {0};

  // Per thread accumulators
  double** txW  = new double*[numThreads];
  double** txD2 = new double*[numThreads];
  long int** tnxc = new long int*[numThreads];
  double** txc = new double*[numThreads];

  double *texpNData = new double[numThreads]{};
  double *teffNData = new double[numThreads]{};

  // Allocate each thread's array
  for (int z = 0; z < numThreads; z++) {
    txW[z] = new double[MAXBINS]{};
    txD2[z] = new double[MAXBINS]{};
    tnxc[z] = new long int[MAXBINS]{};
    txc[z] = new double[MAXBINS]{};
  }

  int* binMaxes = new int[numThreads]{};
  double expNData = 0;
  double effNData = 0;
  int increbinBase;

  if (sampleInfo->hayrecentbins){
    increbinBase=3;
  }
  else{
    increbinBase=0;
  }

  omp_set_num_threads(numThreads);
  double frecI, frecJ;
  std::sort(valid_idx, valid_idx + sampleInfo->numSegLoci);

  std::vector<int> validCromo(sampleInfo->numSegLoci);
  std::vector<int> validChrEnd(sampleInfo->numSegLoci);
  std::vector<double> validPosiCM(sampleInfo->numSegLoci);
  std::vector<uint8_t> chrIsMonotonic(popInfo->numCromo + 1, 1);
  std::vector<uint8_t> validChrIsMonotonic(sampleInfo->numSegLoci);
  for (int locus = 1; locus < popInfo->numLoci; ++locus) {
    const int cromo = popInfo->cromo[locus];
    if (cromo == popInfo->cromo[locus - 1] &&
        popInfo->posiCM[locus] < popInfo->posiCM[locus - 1]) {
      chrIsMonotonic[cromo] = 0;
    }
  }
  for (int idx = 0; idx < sampleInfo->numSegLoci; ++idx) {
    const int locus = valid_idx[idx];
    validCromo[idx] = popInfo->cromo[locus];
    validPosiCM[idx] = popInfo->posiCM[locus];
  }
  for (int start = 0; start < sampleInfo->numSegLoci;) {
    int end = start + 1;
    while (end < sampleInfo->numSegLoci &&
           validCromo[end] == validCromo[start]) {
      ++end;
    }
    const uint8_t monotonic = chrIsMonotonic[validCromo[start]];
    for (int idx = start; idx < end; ++idx) {
      validChrEnd[idx] = end;
      validChrIsMonotonic[idx] = monotonic;
    }
    start = end;
  }

  sampleInfo->numSNPs = 0;
  int repes=1;
  if (sampleInfo->haplotype==3){
    repes=10;
  }
  const bool useBitsetD2 = (repes == 1);
  const int bitWords = (sampleInfo->sampleSizeIndi + 63) / 64;
  // Three packed bitsets per segregating locus:
  //   hetBits    — sample bit set iff genotype == 1
  //   homAltBits — sample bit set iff genotype == 2
  //   calledBits — sample bit set iff genotype < 9 (i.e. not missing)
  // calledBits is the positive of the missing-mask AND carries the
  // sample-size mask baked in: the encode loop only iterates up to
  // sampleSizeIndi, so bits beyond the sample range stay zero in
  // every bitset and no lastWordMask is needed in the inner loop.
  std::vector<uint64_t> hetBits;
  std::vector<uint64_t> homAltBits;
  std::vector<uint64_t> calledBits;
  if (useBitsetD2) {
    const size_t bitsetSize =
        static_cast<size_t>(sampleInfo->numSegLoci) * bitWords;
    hetBits.assign(bitsetSize, 0);
    homAltBits.assign(bitsetSize, 0);
    calledBits.assign(bitsetSize, 0);
    for (int idx = 0; idx < sampleInfo->numSegLoci; ++idx) {
      const int locus = valid_idx[idx];
      const size_t baseOffset = static_cast<size_t>(idx) * bitWords;
      for (int sampleIdx = 0; sampleIdx < sampleInfo->sampleSizeIndi;
           ++sampleIdx) {
        const int indi = sampleInfo->shuffledIndi[sampleIdx];
        const uint8_t genotype = popInfo->indi[indi][locus];
        const size_t wordOffset = baseOffset + sampleIdx / 64;
        const uint64_t bit = uint64_t{1} << (sampleIdx & 63);
        if (genotype == 1) {
          hetBits[wordOffset] |= bit;
          calledBits[wordOffset] |= bit;
        } else if (genotype == 2) {
          homAltBits[wordOffset] |= bit;
          calledBits[wordOffset] |= bit;
        } else if (genotype < 9) {
          // genotype 0 — homozygous reference, called but neither het nor altHom.
          calledBits[wordOffset] |= bit;
        }
        // genotype ≥ 9 → missing; leave all three bits clear.
      }
    }
  }
  const double maxUsefulCenmor = sampleInfo->hayrecentbins ? 56.0 : supx;
  long int completedLoci = 0;
  #pragma omp parallel for private(D, W, D2, W2, recrate, generacion, bin, extrabin,i, j, _i, genot1, genot2, locus1, locus2, tacuHoHo, tacuHoHetHetHo, tacuHetHet, pIndi, cenmor, contaIndX, frecI, frecJ, ss,_irepe) schedule(dynamic, 1000)
  for (i = 0; i < sampleInfo->numSegLoci - 1; ++i) {
    int tid = omp_get_thread_num();
    locus1 = valid_idx[i];
    contaIndX = 0;
    double threadEffNCounter = 0;
    double threadExpNCounter = 0;
    double posiCM1 = validPosiCM[i];
    const int chrEnd = validChrEnd[i];
    const bool cMIsMonotonic = validChrIsMonotonic[i] != 0;
    bool cond;
    int firstJ = i + 1;
    if (cMIsMonotonic) {
      const double minUsefulPosiCM = posiCM1 + infx;
      firstJ = static_cast<int>(
          std::lower_bound(validPosiCM.begin() + firstJ,
                           validPosiCM.begin() + chrEnd,
                           minUsefulPosiCM) -
          validPosiCM.begin());
      while (firstJ > i + 1 &&
             validPosiCM[firstJ - 1] - posiCM1 >= infx) {
        --firstJ;
      }
      while (firstJ < chrEnd && validPosiCM[firstJ] - posiCM1 < infx) {
        ++firstJ;
      }
    }
    for (j = firstJ; j < chrEnd; ++j) {
      locus2 = valid_idx[j];
      cenmor = fabs(validPosiCM[j] - posiCM1);
      if (cMIsMonotonic && cenmor > maxUsefulCenmor) {
        break;
      }
      
      if (sampleInfo->hayrecentbins){
        if ((cenmor > 56)) {continue;}
        if ((cenmor > 30) && (cenmor < 39)) {continue;}
        if ((cenmor > 15) && (cenmor < 20)) {continue;}
        if ((cenmor > supx) && (cenmor < 10)) {continue;}
        if ((cenmor < infx)) {continue;}
        extrabin=-1;
        if (cenmor>=39){extrabin=0;}
        else if(cenmor>=20){extrabin=1;}
        else if(cenmor>=10){extrabin=2;}
      }
      else{
        if ((cenmor > supx) || (cenmor < infx)) {continue;}
        extrabin=-1;
      }

      // acumuladores de genotipos
      D2=W2=0;
      if (useBitsetD2) {
        long int contaIndXInt = 0;
        long int frecIInt = 0;
        long int frecJInt = 0;
        long int tacuHoHoInt = 0;
        long int tacuHoHetHetHoInt = 0;
        long int tacuHetHetInt = 0;
        const size_t offset1 = static_cast<size_t>(i) * bitWords;
        const size_t offset2 = static_cast<size_t>(j) * bitWords;
        for (int word = 0; word < bitWords; ++word) {
          const uint64_t het1    = hetBits   [offset1 + word];
          const uint64_t het2    = hetBits   [offset2 + word];
          const uint64_t hom1    = homAltBits[offset1 + word];
          const uint64_t hom2    = homAltBits[offset2 + word];
          const uint64_t called1 = calledBits[offset1 + word];
          const uint64_t called2 = calledBits[offset2 + word];
          const uint64_t calledBoth = called1 & called2;

          contaIndXInt += PopCount64(calledBoth);
          frecIInt += PopCount64(het1 & called2) +
                      2 * PopCount64(hom1 & called2);
          frecJInt += PopCount64(het2 & called1) +
                      2 * PopCount64(hom2 & called1);
          tacuHetHetInt += PopCount64(het1 & het2);
          tacuHoHetHetHoInt += PopCount64((hom1 & het2) | (het1 & hom2));
          tacuHoHoInt += PopCount64(hom1 & hom2);
        }
        threadExpNCounter += sampleInfo->sampleSizeIndi;
        contaIndX = static_cast<double>(contaIndXInt);
        if (contaIndX == 0) {
          continue;
        }
        frecI = static_cast<double>(frecIInt);
        frecJ = static_cast<double>(frecJInt);
        tacuHoHo = static_cast<double>(tacuHoHoInt);
        tacuHoHetHetHo = static_cast<double>(tacuHoHetHetHoInt);
        tacuHetHet = static_cast<double>(tacuHetHetInt);
        frecI /= 2 * contaIndX;
        frecJ /= 2 * contaIndX;
        W = frecI * frecJ;
        if (params->haplotype == 0){
          D = -2.0 * W +
            (2.0 * tacuHoHo + tacuHoHetHetHo + tacuHetHet / 2.0) / contaIndX;
        }
        else{
          D= tacuHoHo / contaIndX - frecI * frecJ;
        }
        D *= D;
        W *= (1.0 - frecI) * (1.0 - frecJ);
        D2 = D;
        W2 = W;
      } else {
        for (_irepe=0;_irepe<repes;++_irepe){
          tacuHoHo = tacuHoHetHetHo = tacuHetHet = frecI = frecJ = 0;
          pIndi = sampleInfo->shuffledIndi;
          contaIndX = 0;
          for (_i = 0; _i < sampleInfo->sampleSizeIndi; ++_i) {
            genot1=popInfo->indi[*pIndi][locus1];
            genot2=popInfo->indi[*pIndi][locus2];
            if (sampleInfo->haplotype==3){  // Alternative for haplotype=3
              if (genot1==1){
                if (xprng.uniforme01() > 0.5){
                  genot1=0;
                }
                else{
                  genot1=2;
                }
              }
              if (genot2==1){
                if (xprng.uniforme01() > 0.5){
                  genot2=0;
                }
                else{
                  genot2=2;
                }
              }
            }
            ss = genot1 + genot2;
            if (_irepe==0){threadExpNCounter++;}
            cond = (ss < 9);
            contaIndX += cond;
            frecI += genot1 * cond;
            frecJ += genot2 * cond;
            tacuHetHet += (ss == 2) * (genot1 == 1);
            tacuHoHetHetHo += (ss == 3);
            tacuHoHo += (ss == 4); // Homo_no_ref and Homo_no_ref
            ++pIndi;
          }
          if (contaIndX == 0) {
            continue;
          }
          frecI /= 2 * contaIndX;
          frecJ /= 2 * contaIndX;
          W = frecI * frecJ;
          if (params->haplotype == 0){
            D = -2.0 * W +
              (2.0 * tacuHoHo + tacuHoHetHetHo + tacuHetHet / 2.0) / contaIndX;
          }
          else{
            D= tacuHoHo / contaIndX - frecI * frecJ;
          }
          D *= D;
          W *= (1.0 - frecI) * (1.0 - frecJ);
          D2+=D;
          W2+=W;
        }
        D2/=repes;
        W2/=repes;
      }
      if (params->distance == 1) {
        recrate = (1.0 - exp(-0.02 * cenmor)) / 2.0;
      } else if (params->distance == 2) {
        recrate = exp(0.04 * cenmor);
        recrate = 0.5 * (recrate - 1) / (recrate + 1);
      } else {
        recrate = cenmor / 100.0;
        if (recrate > 0.5) {
          recrate = 0.5;
        }
      }
      if (extrabin < 0){
        generacion = 0.5 / ( recrate);
        bin = trunc((generacion - gBase) / 5);
        bin = bin+increbinBase;
      }
      else{
        bin=extrabin;
      }
      // Maybe store binMax in a per-thread array and then loop in the end
      // through all the local binMaxes
      if (bin > binMaxes[tid]) {
        binMaxes[tid] = bin;
      }
      threadEffNCounter += contaIndX;
      //teffNData[tid] += contaIndX;
      tnxc[tid][bin] += 1;
      txc[tid][bin] += 1/recrate; // Para calcular media harmonica
      txD2[tid][bin] += D2;
      txW[tid][bin] += W2;
    }
    teffNData[tid] += threadEffNCounter;
    texpNData[tid] += threadExpNCounter;
    long int progressValue = 0;
#pragma omp atomic capture
    progressValue = ++completedLoci;
    if (progressValue == 1 || progressValue % 100 == 0 ||
        progressValue == sampleInfo->numSegLoci - 1) {
#pragma omp critical(gone_progress_update)
      {
        params->progress.SetStatusDetail(
            "Processed " + std::to_string(progressValue) + " of " +
            std::to_string(sampleInfo->numSegLoci - 1) +
            " segregating loci");
        params->progress.SetTaskProgress(static_cast<float>(progressValue));
        // Mix mode has no GA to publish per-round Ne snapshots, so the
        // chart would stay blank for the whole d² phase. Push a naive
        // Hayes-style curve from the partial accumulators every 2000
        // loci; the Fst-corrected curve replaces it at the end.
        if (params->mix &&
            (progressValue % 2000 == 0 ||
             progressValue == sampleInfo->numSegLoci - 1)) {
          PushNaiveNeFromAccumulators(numThreads, binMaxes,
                                      tnxc, txc, txD2, txW,
                                      params->progress, "Ne_metapop");
        }
      }
    }
  }
  // Reduction of data
  for (int nThread = 0; nThread < numThreads; ++nThread) {
    for (int nBin = 0; nBin <= binMaxes[nThread]; ++nBin) {
      nxc[nBin] += tnxc[nThread][nBin];
      xc[nBin]  += txc[nThread][nBin];
      xD2[nBin] += txD2[nThread][nBin];
      xW[nBin]  += txW[nThread][nBin];
    }

    if (binMaxes[nThread] > sampleInfo->binMax) {
      sampleInfo->binMax = binMaxes[nThread];
    }
    effNData += teffNData[nThread];
    expNData += texpNData[nThread];
  }

  sampleInfo->binExtra= increbinBase;
  refextra=nxc[increbinBase]; // Elimina los bines extras con pocos datos

  for (i=increbinBase-1; i>=0; --i){
      if (nxc[i]<refextra*(3-i) * 0.6){
          for (j=i;j<sampleInfo->binMax;++j){
              nxc[j] = nxc[j+1];
              xc[j]  = xc[j+1];
              xD2[j] = xD2[j+1];
              xW[j]  = xW[j+1];
          }
        --sampleInfo->binMax;
        --sampleInfo->binExtra;
      }
  }

  popInfo->expNData = expNData;
  popInfo->effNData = effNData;
  double XXW=0;
  long int XXn=0;
  for (j = 0; j <= sampleInfo->binMax; ++j) {
    if (nxc[j] > 0) {
      sampleInfo->nxc[j] = nxc[j];
      sampleInfo->xc[j] = xc[j];
      sampleInfo->xc[j] /= sampleInfo->nxc[j];
      sampleInfo->xc[j] = 1 / sampleInfo->xc[j]; // media harmonica
      sampleInfo->numSNPs += sampleInfo->nxc[j];  // numero de parejas de SNPs
      sampleInfo->d2[j] = xD2[j] / xW[j];
        XXW+=xW[j];
        XXn+=nxc[j];
    }
  }

  popInfo->propMiss = 1.0 - static_cast<double>(popInfo->effNData) /
                                static_cast<double>(popInfo->expNData);
  // Numero efectivo de individuos
  sampleInfo->numIndiEff = sampleInfo->sampleSizeIndi * popInfo->effNData / popInfo->expNData;

  sampleInfo->increW=0;
  if (params->basecallerror>0){
    XXW/=XXn;
    epsilon = params->basecallerror;
    iota = (1-2/(2* sampleInfo->numIndiEff +1));
    kappa = 1/(sampleInfo->numIndiEff +1);
    sampleInfo->increW=XXW-pow(((pow(XXW,.5)- epsilon)/(1-4 * epsilon)),2);
    for (j = 0; j <= sampleInfo->binMax; ++j) {
      if (nxc[j] > 0) {
        sampleInfo->d2[j] = (xD2[j]/nxc[j]-sampleInfo->increW*kappa) / \
                            (xW[j]/nxc[j]-sampleInfo->increW*iota);
      }
    }
  }

  if (popInfo->haplotype==0){
    popInfo->f =
      (1.0 + sampleInfo->f * (2.0 * popInfo->avgNumIndiAnalyzed - 1.0)) /
      (2.0 * popInfo->avgNumIndiAnalyzed - 1.0 + sampleInfo->f);
  }
  else if(popInfo->haplotype==2){
    popInfo->f =
     (1.0 + sampleInfo->f * (popInfo->avgNumIndiAnalyzed - 1.0)) /
      (popInfo->avgNumIndiAnalyzed - 1.0 + sampleInfo->f);
  }

  // Free the rest of the memory
  delete[] tnxc;
  delete[] txc;
  delete[] txD2;
  delete[] txW;
  delete[] teffNData;
  delete[] texpNData;
  delete[] binMaxes;
  delete[] valid_idx;
}

void CalculateD2ParallelFst(
  std::string fichero, PopulationInfoMix* popInfoMix,
  PopulationInfo* popInfo, SampleInfo* sampleInfo,
  AppParams* params) {
  // Calculo de D2:

  int* valid_idx = new int[popInfo->numLoci];
  int counter = 0;
  int i, j;
  int _i, _irepe;
  int bin;
  double desde;
  int genot1,genot2;
  bool okchromsize;
  double Dw2, Db2, DbDw;
  double Fst;
  // Initialised here so older gcc's flow analysis doesn't warn —
  // Nemejor is only ever read when Fstmejor < 1 (i.e. an inner
  // iteration actually wrote to it).
  double mind2dif = 0, minm = 0, Nemejor = 0, minmmejor = 0;
  double dw2,db2,dbdw;
  double d2spred1, d2spred2;
  int nsubs;
  bool flag_noconverge=false;

  params->progress.InitCurrentTask(sampleInfo->numSegLoci - 1);

  // Inicializamos la tabla de índices válidos (en los que segrega[x] es true)
  for (int idx = 0; idx < popInfo->numLoci; idx++) {
    if (sampleInfo->segrega[sampleInfo->shuffledLoci[idx]]) {
      valid_idx[counter] = sampleInfo->shuffledLoci[idx];
      counter++;
    }
  }

  if (params->cMMb > 0) {
    for (i = 0; i < popInfo->numLoci; ++i) {
      popInfo->posiCM[i] = params->cMMb * popInfo->posiBP[i] / 1000000;
    }
    popInfo->Mtot = params->cMMb * popInfo->Mbtot / 100;
  }

  // Checking chomosome sizes 
  j=0;
  okchromsize=true;
  desde=popInfo->posiCM[j];
  for (i = 1; i < popInfo->numLoci; ++i) {
    if ((popInfo->cromo[i]!=popInfo->cromo[j]) || (i==(popInfo->numLoci-1))){
      if((popInfo->posiCM[i-1]-desde) < 20.0){
        std::cerr << "Chromosome "<< popInfo->cromo[j]<<" is too small. All chromosomes must be larger than 20cM" << std::endl;
        okchromsize=false;        
      }
      if((popInfo->posiCM[i-1]-desde) > 1000.0){
        std::cerr << "Chromosome "<< popInfo->cromo[j]<<" is too large. All chromosomes must be smaller than 1000cM" << std::endl;
        okchromsize=false;        
      }
      j=i;
      desde=popInfo->posiCM[j];
    }
  }
  if (!okchromsize){
    exit(EXIT_FAILURE);
  }

  double mapdist[MAXDIST] = {0};
  double increcM = 0.01;
  
  int locus1;
  int locus2;
  int* pIndi;

  double cenmor = 0;
  double contaIndX = 0;
  double tacuHoHo, tacuHoHetHetHo, tacuHetHet;
  uint8_t ss;

  double D, W,D2,W2;

  int numThreads = params->numThreads;
  if (numThreads == 0) {
    numThreads = omp_get_max_threads();
  }

  double xW[2] = {0}, xD2[2] = {0};
  long int nxc[2] = {0};

  // Per thread accumulators
  double** txW  = new double*[numThreads];
  double** txD2 = new double*[numThreads];
  long int** tnxc = new long int*[numThreads];
  double** txc = new double*[numThreads];

  double *texpNData = new double[numThreads]{};
  double *teffNData = new double[numThreads]{};

  // Allocate each thread's array
  for (int z = 0; z < numThreads; z++) {
    txW[z] = new double[2]{};
    txD2[z] = new double[2]{};
    tnxc[z] = new long int[2]{};
    txc[z] = new double[2]{};
  }

  double expNData = 0;
  double effNData = 0;

  omp_set_num_threads(numThreads);
  double frecI, frecJ;
  std::sort(valid_idx, valid_idx + sampleInfo->numSegLoci);
  // sampleInfo->numSNPs = 0;
  int repes=1;
  long int completedLoci = 0;
  #pragma omp parallel for private(D, W, D2, W2, bin, i, j, _i, genot1, genot2, locus1, locus2, tacuHoHo, tacuHoHetHetHo, tacuHetHet, pIndi, cenmor, contaIndX, frecI, frecJ, ss,_irepe) schedule(dynamic, 1000)
  for (i = 0; i < sampleInfo->numSegLoci - 1; ++i) {
    int tid = omp_get_thread_num();
    locus1 = valid_idx[i];
    contaIndX = 0;
    double threadEffNCounter = 0;
    double threadExpNCounter = 0;
    double posiCM1 = popInfo->posiCM[locus1];
    int cromo1 = popInfo->cromo[locus1];
    bool cond;
    for (j = i + 1; j < sampleInfo->numSegLoci; ++j) {
      locus2 = valid_idx[j];
      // Are we in the same cromosome?
      bin=-1;
      if (popInfo->cromo[locus2] != cromo1) {
        bin=0;
      }
      else{
        cenmor = fabs(popInfo->posiCM[locus2] - posiCM1);
        if (cenmor>5){ // No se consideran las parejas más próximas de 5 cM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          cenmor/=increcM;
          #pragma omp atomic
          ++mapdist[int(cenmor)];
          bin=1;
        }
      }
      if (bin>=0){
        // acumuladores de genotipos
        D2=W2=0;
        for (_irepe=0;_irepe<repes;++_irepe){
          tacuHoHo = tacuHoHetHetHo = tacuHetHet = frecI = frecJ = 0;
          pIndi = sampleInfo->shuffledIndi;
          contaIndX = 0;
          for (_i = 0; _i < sampleInfo->sampleSizeIndi; ++_i) {
            genot1=popInfo->indi[*pIndi][locus1];
            genot2=popInfo->indi[*pIndi][locus2];
            if (sampleInfo->haplotype==3){  // Alternative for haplotype=3
              if (genot1==1){
                if (xprng.uniforme01() > 0.5){
                  genot1=0;
                }
                else{
                  genot1=2;
                }
              }
              if (genot2==1){
                if (xprng.uniforme01() > 0.5){
                  genot2=0;
                }
                else{
                  genot2=2;
                }
              }
            }
            ss = genot1 + genot2;
            if (_irepe==0){threadExpNCounter++;}
            //++texpNData[tid];
            cond = (ss < 9);
            contaIndX += cond;
            frecI += genot1 * cond;
            frecJ += genot2 * cond;
            tacuHetHet += (ss == 2) * (genot1 == 1);
            tacuHoHetHetHo += (ss == 3);
            tacuHoHo += (ss == 4); // Homo_no_ref and Homo_no_ref
            ++pIndi;
          }
          if (contaIndX == 0) {
            continue;
          }
          frecI /= 2 * contaIndX;
          frecJ /= 2 * contaIndX;
          W = frecI * frecJ;
          if (params->haplotype == 0){
            D = -2.0 * W +
              (2.0 * tacuHoHo + tacuHoHetHetHo + tacuHetHet / 2.0) / contaIndX;
          }
          else{
            D= tacuHoHo / contaIndX - frecI * frecJ;
          }
          D *= D;
          W *= (1.0 - frecI) * (1.0 - frecJ);
          D2+=D;
          W2+=W;
        }
        D2/=repes;
        W2/=repes;

        threadEffNCounter += contaIndX;
        tnxc[tid][bin] += 1;
        txD2[tid][bin] += D2;
        txW[tid][bin] += W2;
      }
    }
    teffNData[tid] += threadEffNCounter;
    texpNData[tid] += threadExpNCounter;
    long int progressValue = 0;
#pragma omp atomic capture
    progressValue = ++completedLoci;
    if (progressValue == 1 || progressValue % 100 == 0 ||
        progressValue == sampleInfo->numSegLoci - 1) {
#pragma omp critical(gone_progress_update)
      {
        params->progress.SetStatusDetail(
            "Processed " + std::to_string(progressValue) + " of " +
            std::to_string(sampleInfo->numSegLoci - 1) +
            " segregating loci");
        params->progress.SetTaskProgress(static_cast<float>(progressValue));
        // Every 200 loci, refresh the corner score box with a partial
        // Fst estimate from the cross-chrom d² accumulator. The Ne
        // chart itself can't refine during this loop — sampleInfo->d2
        // is frozen from CalculateD2Parallel and only the Fst factor
        // moves, by a fraction of a percent on the log-scale chart —
        // so we route the live signal into the score box where it's
        // actually visible.
        if (progressValue % 200 == 0 ||
            progressValue == sampleInfo->numSegLoci - 1) {
          double xD2_part = 0, xW_part = 0;
          long int n_part = 0;
          for (int t = 0; t < numThreads; ++t) {
            xD2_part += txD2[t][0];
            xW_part  += txW[t][0];
            n_part   += tnxc[t][0];
          }
          if (n_part > 100 && xW_part > 0) {
            const double d2s05_partial = xD2_part / xW_part;
            double Fst_partial =
                (d2s05_partial > 0) ? std::sqrt(d2s05_partial) : 0.0;
            if (Fst_partial > 0.9) Fst_partial = 0.9;
            const long int totalLoci = sampleInfo->numSegLoci - 1;
            const int pct = (totalLoci > 0)
                ? static_cast<int>(100.0 * progressValue / totalLoci)
                : 100;
            params->progress.SetBestScoreLabeled(
                Fst_partial, "Best Fst",
                "convergence " + std::to_string(pct) + "%");
          }
        }
      }
    }
  }
  // Reduction of data
  for (int nThread = 0; nThread < numThreads; ++nThread) {
    for (int nBin = 0; nBin <= 1; ++nBin) {
      nxc[nBin] += tnxc[nThread][nBin];
      xD2[nBin] += txD2[nThread][nBin];
      xW[nBin]  += txW[nThread][nBin];
    }
    effNData += teffNData[nThread];
    expNData += texpNData[nThread];
  }
  double numIndi = sampleInfo->sampleSizeIndi * effNData / expNData;
  double effeneind = numIndi;
  
  double xWt = (xW[0]+xW[1])/(nxc[0]+nxc[1]);
  for (j = 0; j <= 1; ++j) {
    if (nxc[j] > 0) {
      xW[j] /= nxc[j];
      xD2[j] /= nxc[j];
    }
  }

  // Calculo de la distribución de distancias para la integral:

  double contadist=0;
  double maxdistance=0;
  int maxdistanceindx=0;
  for (j = 0; j < MAXDIST; ++j){
    if (mapdist[j]>0){
      contadist+=mapdist[j];
      if(mapdist[j]>maxdistance){
        maxdistance=float(j)/10000; // en Morgans
        maxdistanceindx=j;
      }
    }
  }
  for (j = 0; j <=maxdistanceindx; ++j){
      mapdist[j]/=contadist;
  }

  if (popInfo->haplotype==0){
    popInfo->f =
      (1.0 + sampleInfo->f * (2.0 * numIndi - 1.0)) /
      (2.0 * numIndi - 1.0 + sampleInfo->f);
      xWt /= (1-2/(2*numIndi+1));
  }
  else if(popInfo->haplotype==2){
    popInfo->f =
     (1.0 + sampleInfo->f * (numIndi - 1.0)) /
      (numIndi - 1.0 + sampleInfo->f);
      xWt /= (1-2/(numIndi+1));
  }


  double d2s05 = xD2[0]/xW[0];
  double d2slink = xD2[1]/xW[1];
  double fp =popInfo->f;

  double fp12 = (1 + fp) * (1 + fp);
  double sample_size_h = numIndi * 2;
  double samplex = (sample_size_h - 2.0) * (sample_size_h - 2.0) * (sample_size_h - 2.0);
  samplex += 8.0 / 5.0 * (sample_size_h - 2.0) * (sample_size_h - 2.0);
  samplex += 4 * (sample_size_h - 2.0);
  samplex /= ((sample_size_h - 1.0) * (sample_size_h - 1.0) * (sample_size_h - 1.0) + (sample_size_h - 1.0) * (sample_size_h - 1.0));
  double sampley = (2.0 * sample_size_h - 4.0) / ((sample_size_h - 1.0) * (sample_size_h - 1.0));

  double d2 = 0;
  switch (sampleInfo->haplotype){
    case 0:
      d2= (d2s05 - sampley * fp12) / samplex;
      break;
    case 1:
      d2=(d2s05 - 1/(sample_size_h-1));
      break;
    case 2:
      d2=(d2s05 - fp12/(sample_size_h-1));
      break;
    case 3:
      d2=(d2s05 - 1/(numIndi-1));
      break;
  }

  if ((d2>0) && (d2slink>d2s05)){
    // first, find the maximum possible value of Fst
    int ps=2;
    double m;
    double Ne, NeT;
    double KK=0;
    double increFst;
    double maxFst;
    double minFst;
    double d2p05;
    double m2, m3;
    double Fst2, Fst3;
    double Ne2, Ne3;
    double ps2, ps3;
    double minNedif;
    double d2b05= (d2s05 - sampley * fp12) / samplex / 4; // el 4 es por ser unphased
    double kk1=float(ps)/(float(ps)-1);
    double kk2=kk1*kk1;
    double kk4=kk2*kk2;

    if (d2b05>0){
      for (int repe=0;repe<3;++repe){
        int CIEN=100;
        minFst=0;
        maxFst=std::sqrt(d2b05*(ps-1));
        if (maxFst>0.9){maxFst=0.9;}
        increFst=(maxFst-minFst)/CIEN;
        double Fstmejor=999999999,difmejor=999999999;
        for (int ii=0;ii<2;++ii){
          Fst=minFst;
          while (Fst<maxFst){
            Fst+=increFst;
            Mix05ylink(&ps,&kk2,&kk4,&d2slink,&d2s05,&mind2dif,
              &minNedif,&minm, &dw2, &db2, &dbdw,
              &d2spred2,&d2spred1,maxdistance,
              effeneind, fp, 
              Fst, &m, 
              &Dw2,&Db2, &DbDw,&mapdist[0], &Ne,increFst);

            if((Ne<999999999) && (mind2dif<difmejor) && (minm<0.45) && (minm>0)){
              difmejor=mind2dif;
              Nemejor=Ne;
              minmmejor=minm;
              Fstmejor=Fst;
            }
          }
          if (Fstmejor<1){
            minFst=Fstmejor-increFst;
            maxFst=Fstmejor+increFst;
            increFst=(maxFst-minFst)/CIEN;
          }
          else{
            break;
          }
        }
        if (Fstmejor<1){ // Si Fst valido
          Fst=Fstmejor;
          Ne=(Nemejor);
          NeT=Ne*ps;
          // FORMULA ABREVIADA:
          m= (1-Fst)/(4*Ne*Fst*pow((float(ps)/float(ps-1)),2));
          // FORMULA COMPLETA:
          // m = (1-pow(1-1/(2*ps/(ps-1)*Ne)*(1-Fst)/Fst,0.5))/(ps/(ps-1));
          m2=m;
          Fst2=Fst;
          Ne2=Ne;
          ps2=ps;

          double paramssanterior=ps;
          if (repe<2){ // Nueva busqueda si es la primera o segunda pasada
            CIEN=100;
            d2p05=(d2s05-sampley*fp12)/samplex;
            double Fis,lowestFis=999999999;
            bool nopase=true;

            for (nsubs=2;nsubs<101;++nsubs){
              Ne=NeT/nsubs;
              if (Ne>5){
                Fstmejor=999999999;
                double mind2difant=999999999;
                minFst=0;
                maxFst=std::sqrt(d2b05*(nsubs-1));
                if (maxFst>0.9){maxFst=0.9;}
                increFst=(maxFst-minFst)/CIEN;
                for (int ii=0;ii<2;++ii){
                  Fst=minFst;
                  while (Fst<maxFst){
                      Fst+=increFst;
                      CalculaOtros(&mind2dif,&m, Fst, Ne, &nsubs,
                        &Dw2,&Db2, &DbDw, &d2spred2, &ps, &kk2,&kk4, &d2p05);
                      if (mind2dif<1){
                        if (mind2dif<mind2difant){
                          minmmejor=m;
                          mind2difant=mind2dif;
                          Fstmejor=Fst;
                          // mejord2spred2=d2spred2;
                          // mejord2spred1=d2spred1;
                        }
                      }
                    }
                  if (Fstmejor<1){
                    minFst=Fstmejor-increFst;
                    maxFst=Fstmejor+increFst;
                    increFst=(maxFst-minFst)/CIEN;
                  }
                  else{
                    break;
                  }
                }
                Fis=1-(1-fp)/(1-Fstmejor);
                if ((Fstmejor<1) && (Fstmejor>0) && (minmmejor>0) && (minmmejor<0.5)){
                  if (lowestFis<std::abs(Fis)){
                    if (nopase){
                      nopase=false;
                      ps=nsubs-1;
                      break;
                    }
                  }
                  else{
                    lowestFis=std::abs(Fis);
                  }
                }
              }
            }

            if(paramssanterior==ps){ // Stop si la salida de la nueva busqueda no cambia.
              break;
            }
          }
        }
        else{
          flag_noconverge=true;
        }
        if (repe==0){ // Si no converge a la primera, sale
          if (flag_noconverge){
          std::cerr  << "  There is no convergence to a solution within the range of Fst values from 0 to 1.\n";
          std::cerr << "\n";
          exit(EXIT_FAILURE);
          }
        }
        if (repe<2){
          if (!flag_noconverge){ // guarda salida en el buffer salida3 si converge
            m3=m2;
            Fst3=Fst2;
            Ne3=Ne2;
            ps3=ps2;
          }
        }
        if ((ps==2) || (flag_noconverge)){ // recupera buffer anterior y sale
          m2=m3;
          Fst2=Fst3;
          Ne2=Ne3;
          ps2=ps3;
          break;
        }
      }

      m=m2;
      Fst=Fst2;
      Ne=Ne2;
      ps=ps2;

      popInfo->m=m;
      popInfo->Fst=Fst;

      double d2spred=0;
      Ne = FixMixEcuacion05link(&ps,&kk2,&kk4,numIndi, popInfo->f, popInfo->haplotype, popInfo->basecallcorrec, Fst,&popInfo->m,&xWt,&Dw2,&Db2,&DbDw,0.5,&d2spred,d2s05, KK);
      popInfoMix->xc[0]=0.5;
      popInfoMix->Nelink[0]=ps*Ne;
      popInfoMix->mig[0]=popInfo->m;
      popInfoMix->Fst[0]=popInfo->Fst;
      popInfoMix->Dw2[0]=Dw2;
      popInfoMix->Db2[0]=Db2;
      popInfoMix->DbDw[0]=DbDw;
      popInfoMix->Dt2[0]=Dw2+Db2+2*DbDw;
      popInfoMix->Wt[0]=xWt;
      popInfoMix->d2spred[0]=d2spred;
      popInfoMix->d2sobs[0]=d2s05;

      std::string fichero_sal_NeH = fichero + "_GONE2_Ne_mix";
      std::ofstream salida;
      salida.open(fichero_sal_NeH, std::ios::out);
      salida << "# Estimate of historical Ne using observed LD  of different recombination bis. This is an\n"
        "# extension to a subdivided population model of the Hayes et al. (2003) approach and\n"
        "# should be considered an approximate trend.\n"
        "#\n"
        "# Number of subpopulations:\n"
        << ps << "\n"
        "# Fst:\n"
        << std::fixed << std::setprecision(8) << Fst << "\n"
        "# Migration rate:\n"
        << m << "\n"
        "#\n";

      salida <<"# Rec_rate_bin\tgeneration\tN_T_metapop\tNe_metapop\td²_metapop\n";

      // Accumulate Ne_metapop by generation for the live TUI chart.
      // -x mode doesn't run the GA so the chart never receives the
      // per-round SetNeSnapshot pushes — without this it would stay
      // on "waiting for the first GA round…" for the whole run.
      constexpr int kChartMaxGen = 150;
      std::vector<double> chartNe(kChartMaxGen, 0.0);

      double d2sample=d2s05;
      double c=0.5;
      for (int nBin=-1;nBin<sampleInfo->binMax;++nBin){
        if (nBin>=0){
          c=sampleInfo->xc[nBin];
          d2sample=sampleInfo->d2[nBin];
        }
        double generacion=1/(2*c);
        Fst=popInfo->Fst;
        double zDb2 = Fst * Fst / (ps-1);
        m= popInfo->m;
        double ms = float(ps) * m / float(ps-1);
        double m12 = (1-ms)*(1-ms);
        double c2 = c*c;
        double c12 = (1-c) * (1-c);

        double zDw2 = (1-2*Fst*(1-Fst)-zDb2)*(1+c2)/(2*ps*Ne*(1-c12)+2.2*c12);
        double zDbDw = ((kk2)*(zDb2)*m12*(1-(m))*((m))+(zDw2)/(2*ps*Ne))/(1-m12*(1-c)*(1-1.0/ps/Ne));

        //double zDw2p = (zDw2)*c12 + (kk4)*(zDb2)*(1-(m))*(1-(m))*(m)*(m) + 2*(kk2)*(zDbDw)*(1-(m))*(m)*(1-c);
        double zDbDwp = (zDbDw)*m12*(1-c) + (kk2)*(zDb2)*m12*(1-(m))* (m);
        double zDb2p = (zDb2)*m12*m12;   

        double zDw2preconstruida = (d2sample - sampley*fp12)/samplex-4*(zDb2p+zDbDwp);
        double zDw2reconstruida = (zDw2preconstruida - (kk4)*(zDb2)*(1-(m))*(1-(m))*(m)*(m) - 2*(kk2)*(zDbDw)*(1-(m))*(m)*(1-c)) /c12;
        double d2pob=zDw2reconstruida / (1-2*Fst*(1-Fst)-(zDb2)); // ahora es la poblacional
        Ne=(1+c2-2.2*d2pob*c12)/(2*d2pob*(1-c12));
        if (Ne > 0) {
          salida <<std::fixed<< std::setprecision(5)<< c<<"\t"<<std::fixed<< std::setprecision(0)<< generacion<<"\t"<<std::fixed<< std::setprecision(2)<< Ne<< "\t"<< Ne/(1-Fst)<< "\t";
          salida <<std::fixed<< std::setprecision(8)<< d2pob<<"\n";
          const int genIdx =
              static_cast<int>(std::round(generacion)) - 1;
          if (genIdx >= 0 && genIdx < kChartMaxGen) {
            chartNe[genIdx] = Ne / (1 - Fst);
          }
        }


        #ifdef EXPERIMENTAL_MIX
        if (nBin>=0){
            sampleInfo->d2[nBin]=d2pob;
          }
        #endif
      }
      salida.close();

      // Forward-fill: bins map to specific generations, so adjacent
      // gens may be empty. Step-fill from the previous filled value
      // gives a continuous curve for the chart.
      double lastNe = 0.0;
      for (int g = 0; g < kChartMaxGen; ++g) {
        if (chartNe[g] > 0) lastNe = chartNe[g];
        else if (lastNe > 0) chartNe[g] = lastNe;
      }
      if (lastNe > 0) {
        params->progress.SetNeSnapshot(chartNe, "Ne_metapop");
      }

      #ifdef EXPERIMENTAL_MIX
      int firstLowBin = 0;
      for (int nBin = 0; nBin < sampleInfo->binMax; ++nBin) {
        if (sampleInfo->nxc[nBin] < 0.05) {
          firstLowBin = nBin;
          break;
        }
      }
      for (int nBin = firstLowBin; nBin < sampleInfo->binMax; ++nBin) {
        sampleInfo->nxc[nBin-firstLowBin] = sampleInfo->nxc[nBin];
        sampleInfo->d2[nBin-firstLowBin]  = sampleInfo->d2[nBin];
        sampleInfo->xc[nBin-firstLowBin]  = sampleInfo->xc[nBin];
      }
      sampleInfo->binMax -= firstLowBin;
      #endif
    }
  }
  else{
    if (d2 < 0) {
      std::cerr << "# Ne cannot be estimated. d2 estimate for the population is negative.\n" << std::endl;
                 
    } else {
      std::cerr << "# Ne cannot be estimated because the LD is smaller for linked than for unlinked markers.\n"<<std::endl;
    }
  }
  // Free the rest of the memory
  delete[] tnxc;
  delete[] txc;
  delete[] txD2;
  delete[] txW;
  delete[] teffNData;
  delete[] texpNData;
  delete[] valid_idx;
}


void Mix05ylink(int* ps,double* kk2,double* kk4,double* d2slink, double* d2s05,double* mind2dif,
  double* minNedif, double* minm, double* dw2, double* db2, double* dbdw,
  double* d2spred2,double* d2spred1,double maxdistance,
  double effeneind, double fp, 
  double Fst, double* m, 
  double* Dw2,double* Db2, double* DbDw,double* mapdist, double* Ne,double increFst){

  int ii, j,minj;
  double increL, sumad2spred, sumafrec;
  double d2ppred1,d2ppred2;
  double distmin;
  double minlog, maxlog, increlog,logNe;
  double dif;
  double ms, m12, distancia;
  double Dw2p, Db2p, DbDwp;
  double c, c2, c12;
  double fp12, sample_size, sample_size_h, samplex, sampley;

  fp12 = (1 + fp) * (1 + fp);
  sample_size = effeneind;
  sample_size_h = sample_size * 2;
  samplex = (sample_size_h - 2.0) * (sample_size_h - 2.0) * (sample_size_h - 2.0);
  samplex += 8.0 / 5.0 * (sample_size_h - 2.0) * (sample_size_h - 2.0);
  samplex += 4 * (sample_size_h - 2.0);
  samplex /= ((sample_size_h - 1.0) * (sample_size_h - 1.0) * (sample_size_h - 1.0) + (sample_size_h - 1.0) * (sample_size_h - 1.0));
  sampley = (2.0 * sample_size_h - 4.0) / ((sample_size_h - 1.0) * (sample_size_h - 1.0));

  int CIEN=50;

  increL=0.0001; // El cromosoma de 10 Morgans (max) se parte en incres de 1/100.000 Morgans
  distmin=0.05+increL/2;
  minj=int(distmin/increL);

  minlog=0;
  maxlog=6; //maximum Ne = 10^6 per subpop
  increlog=maxlog/CIEN;
  *mind2dif=99999999; // no valido
  *minNedif=99999999;
  *minm=99999999;
  for (ii=0;ii<3;++ii){
      logNe=minlog;
      while(logNe<maxlog){
          logNe+=increlog;
          *Ne=pow(10,logNe);
          // FORMULA ABREVIADA:
          *m= (1-Fst)/(4*(*Ne)*Fst*pow((float(*ps)/float(*ps-1)),2));
          // FORMULA COMPLETA:
          // m = (1-pow(1-1/(2*ps/(ps-1)*Ne)*(1-Fst)/Fst,0.5))/(ps/(ps-1));
          ms=float(*ps)*(*m)/float(*ps-1);
          m12=(1-ms)*(1-ms);
          if ((*m>0) && (*m<0.5)){
              // AHORA VAN LOS DEl MISMO CROMOSOMA
              sumad2spred = 0;
              // sumad2ppred = 0;
              distancia = distmin;
              sumafrec = 0;
              j = minj;
              while (distancia < maxdistance)
              {
                  c = (1 - exp(-2 * distancia)) / 2;
                  c2 = c * c;
                  c12 = (1 - c) * (1 - c);

                  // FORMULACION ABREVIADA:
                  *Db2=Fst*Fst/((*ps)-1);
                  // *Dw2=(1-Fst)*(1-Fst)*(1+c2)/(2*(*ps)*(*Ne)*(1-c12)+2.2*c12);//****************
                  *Dw2=(1-2*Fst*(1-Fst)-(*Db2))*(1+c2)/(2*(*ps)*(*Ne)*(1-c12)+2.2*c12);
                  // *DbDw=((*kk2)*(*Db2)*(*m))/(1-m12*(1-c));// **************************
                  *DbDw=((*kk2)*(*Db2)*m12*(1-(*m))*(*m)+(*Dw2)/(2*(*ps)*(*Ne)))/(1-m12*(1-c)*(1-1/(*ps)/(*Ne)));

                  Dw2p = (*Dw2)*c12 + (*kk4)*(*Db2)*(1-(*m))*(1-(*m))*(*m)*(*m) + 2*(*kk2)*(*DbDw)*(1-(*m))*(*m)*(1-c);
                  DbDwp = (*DbDw)*m12*(1-c) + (*kk2)*(*Db2)*m12*(1-(*m))* (*m);
                  Db2p = (*Db2)*m12*m12;                
                  d2ppred1=Dw2p+4*Db2p+4*DbDwp;           
                  *d2spred1=(d2ppred1)*samplex+sampley * fp12;

                  sumad2spred += (*(mapdist+j)) * (*d2spred1); // La integral
                  distancia += increL;
                  sumafrec += (*(mapdist+j));
                  ++j;
              }
              if (sumafrec > 0)
              {
                  *d2spred1 = (sumad2spred / sumafrec);
                  // d2ppred1 = (sumad2ppred / sumafrec);
              }
              else
              {
                  std::cerr << "Within chromosome integral cannot be computed." << std::endl;
                  exit(EXIT_FAILURE);
              }
              //
              // AHORA VAN LOS DE DISTINTO CROMOSOMA
              c = 0.5;
              c2 = c * c;
              c12 = (1 - c) * (1 - c);

              // FORMULACION ABREVIADA:
              *Db2=Fst*Fst/((*ps)-1);
              // *Dw2=(1-Fst)*(1-Fst)*(1+c2)/(2*(*ps)*(*Ne)*(1-c12)+2.2*c12);// *********************
              *Dw2=(1-2*Fst*(1-Fst)-(*Db2))*(1+c2)/(2*(*ps)*(*Ne)*(1-c12)+2.2*c12);
              // *DbDw=((*kk2)*(*Db2)*(*m))/(1-m12*(1-c));// **************************
              *DbDw=((*kk2)*(*Db2)*m12*(1-(*m))*((*m))+(*Dw2)/(2*(*ps)*(*Ne)))/(1-m12*(1-c)*(1-1/(*ps)/(*Ne)));

              Dw2p = (*Dw2)*c12 + (*kk4)*(*Db2)*(1-(*m))*(1-(*m))*(*m)*(*m) + 2*(*kk2)*(*DbDw)*(1-(*m))*(*m)*(1-c);
              DbDwp = (*DbDw)*m12*(1-c) + (*kk2)*(*Db2)*m12*(1-(*m))* (*m);
              Db2p = (*Db2)*m12*m12;
              d2ppred2=Dw2p+4*Db2p+4*DbDwp;           
              *d2spred2=(d2ppred2)*samplex+sampley * fp12;

              dif=pow(((*d2spred1)-(*d2slink))/(*d2slink),2)+pow(((*d2spred2)-(*d2s05))/(*d2s05),2);
              if (dif<(*mind2dif)){
                  *mind2dif=dif;
                  *minNedif=(*Ne);
                  *minm=*m;
                  *dw2=*Dw2;
                  *db2=*Db2;
                  *dbdw=*DbDw;
              }
          }
      }
      *Ne=*minNedif;
      if ((*Ne)==99999999){
          break;
      }
      minlog=log10(*Ne)-increlog;
      maxlog=log10(*Ne)+increlog;
      increlog=(maxlog-minlog)/CIEN;    
  }
}

void CalculaOtros(double* mind2dif,double* m, double Fst,double Ne, int* nsubs,
  double* Dw2,double* Db2, double* DbDw, double* d2spred2, int* ps,double* kk2,double* kk4, double* d2p05){

  double Dw2p, Db2p, DbDwp;
  double c, c2, c12;
  double m12, ms;

  c=0.5;
  c2 = c * c;
  c12 = (1 - c) * (1 - c);

  *mind2dif=999999999;
  *m= (1-Fst)/(4*Ne*Fst*pow((float(*nsubs)/float(*nsubs-1)),2));
  ms=float(*nsubs)*(*m)/float(*nsubs-1);
  m12=(1-ms)*(1-ms);


  if (((*m)>0) && (*m<0.5)){

      // FORMULACION ABREVIADA:
      *Db2=Fst*Fst/(*nsubs-1);
      // *Dw2=(1-Fst)*(1-Fst)*(1+c2)/(2*(*nsubs)*Ne*(1-c12)+2.2*c12);// ***************************
      *Dw2=(1-2*Fst*(1-Fst)-(*Db2))*(1+c2)/(2*(*nsubs)*Ne*(1-c12)+2.2*c12);
      // *DbDw=((*kk2)*(*Db2)*(*m))/(1-m12*(1-c));// **************************
      *DbDw=((*kk2)*(*Db2)*m12*(1-(*m))*(*m)+(*Dw2)/(2*(*nsubs)*Ne))/(1-m12*(1-c)*(1-1/(*nsubs)/Ne));

      Dw2p = (*Dw2)*c12 + (*kk4)*(*Db2)*(1-(*m))*(1-(*m))*(*m)*(*m) + 2*(*kk2)*(*DbDw)*(1-(*m))*(*m)*(1-c);
      DbDwp = (*DbDw)*m12*(1-c) + (*kk2)*(*Db2)*m12*(1-(*m))* (*m);
      Db2p = (*Db2)*m12*m12;                
      *d2spred2=(Dw2p+4*Db2p+4*DbDwp); // en realidad es poblacional

      *mind2dif=std::abs(*d2spred2-*d2p05); // diferencia entre esperado y real.
  }
}

//  MIGRACION entre subpoblaciones:
double FixMixEcuacion05link(int* ps,double* kk2,double* kk4,double effeneind, double fp, int haplotype, double basecallcorrec, double Fst, double* m, double* Wt,double* Dw2,double* Db2, double* DbDw, double c,double* d2spred,double d2sobs, double KK)
{
    int DOSMIL, i, ii, conta;
    double sample_size, sample_size_h, samplex, sampley;
    double increNe;
    double fp12, ms;
    bool flagbreak = false;
    double c2,c12,m12;
    double Db2p,Dw2p,DbDwp;

    DOSMIL = 2000;
    fp12 = (1 + fp) * (1 + fp);
    sample_size = effeneind;
    sample_size_h = sample_size * 2;
    samplex = (sample_size_h - 2.0) * (sample_size_h - 2.0) * (sample_size_h - 2.0);
    samplex += 8.0 / 5.0 * (sample_size_h - 2.0) * (sample_size_h - 2.0);
    samplex += 4 * (sample_size_h - 2.0);
    samplex /= ((sample_size_h - 1.0) * (sample_size_h - 1.0) * (sample_size_h - 1.0) + (sample_size_h - 1.0) * (sample_size_h - 1.0));
    sampley = (2.0 * sample_size_h - 4.0) / ((sample_size_h - 1.0) * (sample_size_h - 1.0));

    c2 = c * c;
    c12 = (1 - c) * (1 - c);
    ms=float(*ps)*(*m)/float(*ps-1);
    m12=(1-ms)*(1-ms);

    double Ne = 5000;
    for (ii = 0; ii < 2; ++ii)
    {
        increNe = Ne / (20 * (ii + 1));
        for (i = 0; i < 20; ++i)
        {
            for (conta = 0; conta < DOSMIL; ++conta)
            {

                *Db2=Fst*Fst/((*ps)-1);
                // *Dw2=(1-Fst)*(1-Fst)*(1+c2)/(2*(*ps)*Ne*(1-c12)+2.2*c12);
                *Dw2=(1-2*Fst*(1-Fst)-(*Db2))*(1+c2)/(2*(*ps)*Ne*(1-c12)+2.2*c12);
                // *DbDw=((*kk2)*(*Db2)*(*m))/(1-m12*(1-c));// **************************
                *DbDw=((*kk2)*(*Db2)*m12*(1-(*m))*(*m)+(*Dw2)/(2*(*ps)*Ne))/(1-m12*(1-c)*(1-1/(*ps)/Ne));

                Dw2p = (*Dw2)*c12 + (*kk4)*(*Db2)*(1-(*m))*(1-(*m))*(*m)*(*m) + 2*(*kk2)*(*DbDw)*(1-(*m))*(*m)*(1-c);
                DbDwp = (*DbDw)*m12*(1-c) + (*kk2)*(*Db2)*m12*(1-(*m))* (*m);
                Db2p = (*Db2)*m12*m12;                

                *Dw2=Dw2p;
                *DbDw=DbDwp;
                *Db2=Db2p;
                *d2spred=((*Dw2)+4*(*Db2)+4*(*DbDw))*samplex+sampley*fp12;
      
                if (*d2spred > d2sobs)
                {
                    if (increNe < 0)
                    {
                        increNe /= -5;
                        break;
                    }
                }
                if (*d2spred < d2sobs)
                {
                    if (increNe > 0)
                    {
                        increNe /= -5;
                        break;
                    }
                }
                if ((Ne + increNe) < 3)
                {
                    increNe /= 5;
                    break;
                }
                else
                {
                    Ne += increNe;
                }
                if (abs(increNe) < 0.01)
                {
                    break;
                }
                if (Ne > 100000000)
                {
                    flagbreak = true;
                    break;
                }
            }
            if (abs(increNe) < 0.01)
            {
                break;
            }
            if (flagbreak)
            {
                break;
            }
        }
        if (flagbreak)
        {
            break;
        }
    }
    return Ne;
}

/*
 * GONE2 — Genetic Optimization for Ne Estimation
 *
 * Authors: Enrique Santiago, Carlos Köpke
 */

#include <omp.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

#include "lib/population.hpp"
#include "lib/ped_map.hpp"
#include "lib/tped.hpp"
#include "lib/vcf.hpp"
#include "lib/sample.hpp"
#include "lib/libgone.hpp"
#include "lib/ne_ref.hpp"
#include "lib/kinship.hpp"

#ifdef GONE_NCURSES_TUI
#include "lib/ncurses_tui.hpp"
#endif

int main(int argc, char * argv[]) {
  // RNG used only for shuffling loci and individuals at startup.
  // The GA itself uses xprng (Xoshiro256+) seeded with `params.semilla`.
  std::random_device rd;
  std::mt19937 g(rd());

  // Default-construct in place. ProgressStatus owns a std::mutex,
  // so the older `AppParams params = AppParams();` form fails to
  // move-construct on pre-C++17 (g++ ≤ 16 without guaranteed copy
  // elision) and would need `AppParams params{};` anyway.
  AppParams params;
  HandleInput(argc, argv, &params);
  params.progress.ConfigurePhases({
      "Reading input files",
      "Preparing sample",
      "Calculating d2",
      "Estimating Ne",
      "Writing outputs"});
#if defined(GONE_TUI) && !defined(GONE_NCURSES_TUI)
  params.progress.SetTerminalUiEnabled(!params.quiet);
#else
  params.progress.SetTerminalUiEnabled(false);
#endif
#ifdef GONE_NCURSES_TUI
  if (!params.quiet) {
    StartNcursesTui(&params.progress);
  }
#endif

  // -f reference Ne overlay (loaded for any build; only consumed
  // visually by the ncurses TUI).
  if (!params.realNeFile.empty()) {
    std::vector<double> refNe = LoadReferenceNe(params.realNeFile);
    if (!refNe.empty()) {
      params.progress.SetReferenceNe(refNe);
    }
  }

  int j, j3;
  double start, stop;

  std::string fichProgress = params.fileOut + "_GONE_progress.tmp";

  if (params.numThreads > 0) {
    omp_set_num_threads(params.numThreads);
  }

#if !defined(GONE_TUI) && !defined(GONE_NCURSES_TUI)
  if (!params.quiet) {
    std::cout << " Progress information stored at " << fichProgress << "\n";
    std::cout << " Check it by issuing 'cat " << fichProgress << "'\n";
    std::cout << " Loading data files" << std::endl;
  }
#endif

  start = omp_get_wtime();
  // LECTURA DE LOS DATOS DE SIMULACION:
  // Progress:Phase 1
  params.progress.InitTotalTasks(5, fichProgress.c_str());
  params.progress.SetCurrentTask(0, "Reading input files");
  params.progress.InitCurrentTask(1);
  params.progress.SetStatusDetail(params.fich);
  params.progress.SaveProgress();

  PopulationInfo* popInfo = new PopulationInfo();
  PopulationInfoMix* popInfoMix = new PopulationInfoMix();
  // Allocate individuals
  SampleInfo* sampleInfo = new SampleInfo();

  sampleInfo->haplotype = params.haplotype; // Chapuza dado el lio de declaraciones de variables
  sampleInfo->mix = params.mix;
  popInfo->haplotype = params.haplotype;
  sampleInfo->hayrecentbins = params.hayrecentbins;
  popInfo->hayrecentbins = params.hayrecentbins;
  popInfo->basecallcorrec = 1;

  if (params.ftype == "tped") {
    if (!ReadTped(params.fich, popInfo)) {
      delete popInfo;
      delete sampleInfo;
      exit(EXIT_FAILURE);
    }
  } else if (params.ftype == "vcf") {
    if (!ReadVcf(params.fich, popInfo)) {
      delete popInfo;
      delete sampleInfo;
      exit(EXIT_FAILURE);
    }
    // If vcf check param cMMb
    if (params.cMMb == 0) {
      std::string fichmap = params.fich.substr(0, params.fich.rfind(".")) + ".map";
      if (!ReadMap(fichmap, popInfo)) {
        std::cerr <<
          "\nNo recombination rate has been specified (option -r) "
          "and no map file present (" << fichmap << ")" << std::endl;
        exit(EXIT_FAILURE);
      }
      if (popInfo->Mtot == 0) {
        std::cerr << "\nFound map file \"" << fichmap <<
          "\" but no genetic map was found"
          << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  } else if (params.ftype == "ped") {
    std::string fichped = params.fich;
    std::string fichmap = params.fich.substr(0, params.fich.rfind(".")) + ".map";
    if (!ReadFile(fichped, fichmap, popInfo)) {
      delete popInfo;
      delete sampleInfo;
      exit(EXIT_FAILURE);
    }
  } else {
    std::cerr << "Unknown file extension." << std::endl;
    exit(EXIT_FAILURE);
  }

  params.progress.SetTaskProgress(1);
  double tProcessFile = omp_get_wtime() - start;
#if defined(GONE_TUI) || defined(GONE_NCURSES_TUI)
  (void)tProcessFile;
#endif
#if !defined(GONE_TUI) && !defined(GONE_NCURSES_TUI)
  if (!params.quiet) {
    std::cout << " Reading the input file/s took "
              << std::fixed << std::setprecision(2) << tProcessFile << " sec"
              << std::endl;
  }
#endif

  params.progress.SetCurrentTask(1, "Preparing sample");
  params.progress.InitCurrentTask(4);

  // -k: drop individuals more related than the user-chosen KING-robust
  // φ threshold before any sample stats are computed. No-op when the
  // threshold is 0 (default). Runs early so the surviving sample sees
  // a clean shuffle and downstream allele-frequency / heterozygosity
  // calculations.
  if (params.kinshipThreshold > 0.0) {
    FilterRelatedIndividuals(popInfo, params.kinshipThreshold,
                             &params.progress);
  }

  params.progress.SetStatusDetail("Shuffling loci and individuals");
  params.progress.SaveProgress();

  // Shuffle loci into a random analysis order.
  for (j = 0; j < popInfo->numLoci; ++j) {
    sampleInfo->shuffledLoci[j] = j;
  }
  std::shuffle(sampleInfo->shuffledLoci,
               &(sampleInfo->shuffledLoci[popInfo->numLoci]), g);

  // Shuffle individuals (for haplotype 2 — phased diploids — we shuffle
  // diploid pairs as a unit so the two haplotypes stay adjacent).
  if (popInfo->haplotype != 2) {
    for (j = 0; j < popInfo->numIndi; ++j) {
      sampleInfo->shuffledIndi[j] = j;
    }
    std::shuffle(sampleInfo->shuffledIndi,
                 &(sampleInfo->shuffledIndi[popInfo->numIndi]), g);
  } else {
    const int ndips = popInfo->numIndi / 2;
    for (j = 0; j < ndips; ++j) {
      sampleInfo->shuffledIndi[j] = j;
    }
    std::shuffle(sampleInfo->shuffledIndi,
                 &(sampleInfo->shuffledIndi[ndips]), g);
    // Expand each diploid index into its two haplotype indices,
    // iterating right-to-left to write in place safely.
    for (j = ndips - 1; j >= 0; --j) {
      sampleInfo->shuffledIndi[j * 2]     = sampleInfo->shuffledIndi[j] * 2;
      sampleInfo->shuffledIndi[j * 2 + 1] = sampleInfo->shuffledIndi[j * 2] + 1;
    }
  }

  if (params.numSNPs > 0) {
    sampleInfo->numLoci = params.numSNPs;
  } else {
    sampleInfo->numLoci = popInfo->numLoci;
  }

  if (params.numSample > 0) {
    sampleInfo->sampleSizeIndi = params.numSample;
    if (popInfo->haplotype == 2){sampleInfo->sampleSizeIndi *= 2;}
  } else {
    sampleInfo->sampleSizeIndi = popInfo->numIndi;
  }
  if (sampleInfo->sampleSizeIndi > popInfo->numIndi) {
    sampleInfo->sampleSizeIndi = popInfo->numIndi;
  }
  
  if (sampleInfo->numLoci > popInfo->numLoci) {
    sampleInfo->numLoci = popInfo->numLoci;
  }

  params.progress.SetStatusDetail("Calculating allele frequencies");
  params.progress.SetTaskProgress(1);
  CalculateFrequencies(popInfo, sampleInfo);

  params.progress.SetStatusDetail("Calculating Hardy-Weinberg deviation");
  params.progress.SetTaskProgress(2);
  CalculateF(sampleInfo, popInfo);

  params.progress.SetStatusDetail("Calculating heterozygosity");
  params.progress.SetTaskProgress(3);
  if ((params.haplotype == 0) || (params.haplotype == 2)){
    CalculateHeterozigosity(popInfo, sampleInfo);
  }
  params.progress.SetStatusDetail("Sample preparation complete");
  params.progress.SetTaskProgress(4);

  popInfo->hetEspAll = sampleInfo->hetEspAll;

  // Progress: Phase 4
  params.progress.SetCurrentTask(2, "Calculating d2");
  params.progress.SetStatusDetail("Preparing d2 calculation");
  params.progress.SaveProgress();
#if !defined(GONE_TUI) && !defined(GONE_NCURSES_TUI)
  if (!params.quiet) {
    std::cout << " Measuring d²" << std::endl;
  }
#endif

  CalculateD2Parallel(popInfo, sampleInfo, &params);
 
  std::stringstream salida2;
  salida2
      << params.haplotype
      << "\t# Phase (0:unphased diploids; 1:haploids; 2:phased diploids; 3:low_coverage)\n";
  salida2 << sampleInfo->numIndiEff
          << "\t# Sample size (x2 in phased diploids)\n";
  salida2 << sampleInfo->f << "\t# Hardy-Weinberg deviation in the sample\n";
  salida2 << sampleInfo->binExtra << "\t# Numero de bins extras\n";
  for (j3 = 0; j3 <= sampleInfo->binMax; ++j3) {
    if (sampleInfo->nxc[j3] > 0) {
      salida2 << std::fixed << std::setprecision(0) << sampleInfo->nxc[j3]
              << "\t" << std::fixed << std::setprecision(10)
              << sampleInfo->xc[j3] << "\t" << sampleInfo->d2[j3] << "\n";
    }
  }
  std::string fichsal2 = params.fileOut + "_d2.txt";
  std::ofstream outputFile2;
  outputFile2.open(fichsal2);
  outputFile2 << salida2.str();
  outputFile2.close();

  if (params.cMMb == 0 && popInfo->Mtot == 0) {
    std::cerr << " There is no recombination map." << std::endl;
    delete popInfo;
    delete sampleInfo;
    exit(EXIT_FAILURE);
  }
  // Here cometh GONE!
#if !defined(GONE_TUI) && !defined(GONE_NCURSES_TUI)
  if (!params.quiet) {
    std::cout << " Estimating Ne" << std::endl;
  }
#endif
  // Progress: Phase 5
  params.progress.SetCurrentTask(3, "Estimating Ne");
  params.progress.SetStatusDetail("Preparing GA rounds");
  params.progress.SaveProgress();
  //params.progress.PrintProgress();

  gone(&params, fichsal2, argc, argv, popInfo, sampleInfo);
  if (params.mix) {
    CalculateD2ParallelFst(params.fileOut, popInfoMix, popInfo, sampleInfo,
                           &params);

#ifdef EXPERIMENTAL_MIX
    // Experimental second GA pass on the corrected (population-level)
    // d² values. Disabled by default; rebuild with -DEXPERIMENTAL_MIX
    // to enable.
#if !defined(GONE_TUI) && !defined(GONE_NCURSES_TUI)
    if (!params.quiet) {
      std::cout << " Recalculating Ne" << std::endl;
    }
#endif
    std::remove(fichsal2.c_str());
    std::stringstream salida3;
    salida3 << params.haplotype
            << "\t# Phase (0:unphased diploids; 1:haploids; "
               "2:phased diploids; 3:low_coverage)\n";
    salida3 << sampleInfo->numIndiEff
            << "\t# Sample size (x2 in phased diploids)\n";
    salida3 << sampleInfo->f << "\t# Hardy-Weinberg deviation in the sample\n";
    salida3 << sampleInfo->binExtra << "\t# Number of extra bins\n";
    for (j3 = 0; j3 < sampleInfo->binMax; ++j3) {
      if (sampleInfo->nxc[j3] != 0) {
        salida3 << sampleInfo->nxc[j3] << "\t" << sampleInfo->xc[j3] << "\t"
                << sampleInfo->d2[j3] << "\n";
      }
    }
    std::ofstream outputFile3;
    outputFile3.open(fichsal2);
    outputFile3 << salida3.str();
    outputFile3.close();
    params.mix = false;
    gone(&params, fichsal2, argc, argv, popInfo, sampleInfo);
    params.mix = true;
#endif
  }
  // Remove intermediate d²-stage file.
  std::remove(fichsal2.c_str());

  params.progress.SetCurrentTask(4, "Writing outputs");
  params.progress.InitCurrentTask(2);
  params.progress.SetStatusDetail(
      "Writing statistics and cleaning temporary files");
  params.progress.SaveProgress();

  stop = omp_get_wtime() - start;
  std::stringstream salida;
  salida << "# (GONE v2.0)\n";
  salida << "# Command:";
  for (int i = 0; i < argc; ++i) {
    salida << " " << argv[i];
  }
  salida << "\n";
  salida << "# Running time:";
  salida << static_cast<float>(stop) << "sec\n";
  salida << "#\n";

  salida << "# PREPROCESSING INFORMATION:\n";
  if (params.haplotype == 0){
      salida << "# Temporal Ne estimation using unphased diploid data under the assumption of\n";
  }
  else if (params.haplotype == 1){
      salida << "# Temporal Ne estimation using haploid data under the assumption of\n";
  }
  else if (params.haplotype == 2){
      salida << "# Temporal Ne estimation using phased diploid data under the assumption of\n";
  }
  else if (params.haplotype == 3){
      salida << "# Temporal Ne estimation using diploid low-coverage data under the assumption of\n";
  }
  if (params.mix){
      salida << "#   a metapopulation model with panmixia within subpopulations.\n";
  }
  else{
      salida << "#   a single population model with panmixia.\n";
  }
  if (params.cMMb > 0){
      salida << "# The marker locations in a physical map are known and the chromosome sizes and the\n";
      salida << "#   total genome size will be calculated using those locations. The physical map will be\n";
      salida << "#   converted into a genetic map using a constant recombination rate across the genome: \n";
      salida << "#   distances in Morgans will be calculated using the physical distances (in Mb).\n";
      salida << "# Ne will be inferred using the recombination rates between loci and the weighted cuadratic\n";
      salida << "#   correlations of alleles of loci pairs.\n";
  } 
  else{
      salida << "# The marker locations in a genetic map are known and the chromosome sizes and \n";
      salida << "#   the total genetic size will be calculated using those locations. The genetic \n";
      salida << "#   distances (in Morgans) between loci pairs will be calculated directly using\n";
      salida << "#   the locations in the genetic map.\n";
      salida << "# Ne will be inferred using the recombination rates between loci and the weighted cuadratic\n";
      salida << "#   correlations of alleles of loci pairs.\n";
  }
  if (params.mix){
      salida << "# The historical Ne series will inferred by assuming that observed correlation between\n";
      salida << "#   sites with recombination rate c reflects the Ne of 1/(2c) generations ago.\n#\n";
  }
  else{
      salida << "# The complete theoretical model will be used to search for the historical Ne series\n";
      salida << "#   that best fit the the observed correlations across bins of recomination.\n#\n";
  }
  salida << "# INPUT PARAMETERS:\n";
  salida << "# Type of genotyping data. 0:unphased diploids, 1:haploids, 2:phased diploids, 3:pseudohaploids:\n";
  salida << std::fixed << std::setprecision(0);
  salida << params.haplotype << "\n";
  if ((params.haplotype == 0) || (params.haplotype == 2)){
    if (params.coverage>1){
      salida << "# Coverage (depth of DNA sequencing):\n";
      salida << std::fixed << std::setprecision(1);
      salida << params.coverage << "\n";
    }
  }
  salida << "# Rate of base call errors:\n";
  salida << std::fixed << std::setprecision(5);
  salida << params.basecallerror << "\n";
  salida << "# Number of chromosomes:\n";
  salida << std::fixed << std::setprecision(0);
  salida << popInfo->numCromo << "\n";
  salida << "# Genome size in Morgans:\n";
  salida << std::fixed << std::setprecision(4);
  salida << popInfo->Mtot << "\n";
  salida << "# Genome size in Mb:\n";
  salida << std::fixed << std::setprecision(2);
  salida << popInfo->Mbtot << "\n";
  salida << "# Total number of individuals in the input file:\n";
  salida << std::fixed << std::setprecision(0);
  if (sampleInfo->haplotype == 2){
    salida << popInfo->numIndi/2 << "\n";
  }
  else{
    salida << popInfo->numIndi << "\n";
  }
  salida << "# Number of individuals included in the analysis:\n";
  if (sampleInfo->haplotype == 2){
    salida << sampleInfo->sampleSizeIndi/2 << "\n";
  }
  else{
    salida << sampleInfo->sampleSizeIndi << "\n";
  }
  salida << "# Average Number of individuals included in the analysis (corrected for missing genotypes):\n";
  salida << std::fixed << std::setprecision(4);
  if (sampleInfo->haplotype == 2){
    salida << popInfo->avgNumIndiAnalyzed/2 <<"\n";
  }
  else{
    salida << popInfo->avgNumIndiAnalyzed <<"\n";
  }
  salida << "# Effective Number of individuals for correlations (corrected for missing genotypes):\n";
  if (sampleInfo->haplotype == 2){
    salida << sampleInfo->numIndiEff/2 <<"\n";
  }
  else{
    salida << sampleInfo->numIndiEff <<"\n";
  }
  salida << "# Number of markers in the input file:\n"; //Anadido "markers"
  salida << std::fixed << std::setprecision(0);
  salida << popInfo->numLoci << "\n";
//  salida << "# Number of SNPs included in the analysis:\n"; //Esto sobra
//  salida << std::fixed << std::setprecision(0);
//  salida << sampleInfo->numLoci << "\n";
  salida << "# Number of SNPs included in the analysis (only polymorphic and with less than 20% missing data):\n";
  salida << std::fixed << std::setprecision(0);
  salida << sampleInfo->numSegLoci << "\n";
  salida << std::fixed << std::setprecision(8);
  salida << "# Proportion of missing data:\n";  // DIFIERE DEL preGONE como consecuencia de los dos anteriores
  salida << popInfo->propMiss << "\n";
  salida << "#\n";

  if ((params.haplotype == 0) || (params.haplotype == 2)){
    salida << "# Estimated Fis value of the population (deviation from H-W "
            "proportions):\n";
    salida << popInfo->f << "\n";
    salida << "# Heterozygosity observed in the sample (only polymorphic sites):\n";
    salida << sampleInfo->hetAvg<< "\n";
    salida << "# Heterozygosity observed in the sample (all sites):\n";
    salida << sampleInfo->hetAvgAll<< "\n";
  }
  salida << "# Heterozygosity expected (H-W eq.) in the sample (only polymorphic sites):\n";
  salida << sampleInfo->hetEsp<< "\n";
  salida << "# Heterozygosity expected (H-W eq.) in the sample (all sites):\n";
  salida << sampleInfo->hetEspAll<< "\n";
  params.progress.SetTaskProgress(2);
#ifdef GONE_NCURSES_TUI
  // Write the STATS file *before* tearing down ncurses so the chart
  // stays on screen while the user inspects the result. Hold the
  // screen until Ctrl-C is pressed.
  if (!params.quiet && !params.printToStdOut) {
    std::string fichsal = params.fileOut + "_GONE2_STATS";
    std::ofstream outputFile(fichsal);
    outputFile << salida.str();
  }
  if (!params.quiet) {
    params.progress.SetStatusDetail(
        "Done. Outputs written to " + params.fileOut + "_GONE2_{Ne,d2,STATS}. "
        "Press Ctrl-C to exit.");
    params.progress.SaveProgress();
    WaitForSigint();
  }
  StopNcursesTui();
  if (params.printToStdOut && !params.quiet) {
    std::cout << salida.str() << std::endl;
  }
#else
  if (!params.quiet) {
    if (params.printToStdOut) {
      std::string output = salida.str();
      std::cout << output << std::endl;
    } else {
      std::string fichsal = params.fileOut + "_GONE2_STATS";
      std::ofstream outputFile;
      outputFile.open(fichsal);
      outputFile << salida.str();
      outputFile.close();
    }
  }
#endif

  delete sampleInfo;
  delete popInfo;
  std::remove(fichProgress.c_str());
  return 0;
}

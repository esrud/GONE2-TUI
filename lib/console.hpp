#pragma once

#include <getopt.h>
#include <omp.h>
#include <iostream>
#include <string>
#include <math.h>

#include "constants.hpp"
#include "progress.hpp"

typedef struct {
    int haplotype;
    bool mix;
    double basecallerror;
    double miss;
    double coverage;
    double ngensampling;
    int numThreads;
    int numSample;
    long int numSNPs;
    double hc;
    double lc;
    double cMMb;
    double MAF;
    int distance;
    bool quiet;
    bool printToStdOut;
    std::string fich;
    std::string fileOut;
    std::string ftype;
    // Drop individuals whose pairwise KING-robust kinship exceeds
    // this threshold before the d² calculation runs. 0 = filter off.
    // Common cutoffs: 0.0884 (2nd-degree+), 0.0442 (3rd-degree+).
    double kinshipThreshold;
    // Per-site per-generation mutation rate (-w). When > 0, libgone
    // computes a mutation-drift-balance Ne target from observed
    // heterozygosity and adds a soft penalty on the GA's deep-end
    // plateau so it doesn't drift to physically implausible values.
    // Typical values: ~1e-8 for most vertebrates. 0 = anchor off.
    double mutationRate;
    // Optional path to a reference Ne-history file (from -f).
    // Two accepted formats: simu-param style with a trailing
    // "indiv_diploides generaciones" block, or per-generation
    // "Gen Ne" table (data/simus/data_*.true). The ncurses TUI
    // overlays this curve in red so the user can eyeball the GA
    // result against ground truth.
    std::string realNeFile;
    int flags;
    int muestraSalida;
    int semilla;
    int sizeBins;
    int nbins;
    bool hayrecentbins;
    ProgressStatus progress;
} AppParams;


void HandleInput(int argc, char * argv[], AppParams* params);
void SetDefaultParameters(AppParams* params);
bool GetFileType(std::string fname, std::string *ftype);

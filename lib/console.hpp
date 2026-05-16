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

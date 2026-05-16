#pragma once

#include <algorithm>
#include <string>
#include "./bicho.hpp"
#include "./ga_config.hpp"

typedef struct {
  int gmax[3];
  int topePosInicio[4];
  double fval;
} PoolParams;

typedef struct {
  double frecMut;
  double efectoMut;
  double efectoMutLateral;
  double frecRec;
  double frecNoRec;
  double frecMutLateral;
  int maxSegmentos;
  bool solape;
} MutParams;

typedef struct {
  Bicho parents[MAX_POPULATION];
  Bicho children[MAX_POPULATION];
  PoolParams poolParams;
  MutParams mutParams;
  int numGenerations;
  int numChildren;
  // Runtime knob consulted by MutateNeRnd. Set by Run() per
  // generation under GA_TWO_STAGE, or constant under GA_HEAVY_TAIL.
  bool heavyTailActive;
  // Kick / scout knobs read by MutateNeRnd and Run. Initialised
  // from the compile-time -DGA_* defaults in libgone.cpp.
  GAConfig gaConfig;
} Pool;

void SetInitialPoolParameters(Pool* pool, double mincval);
void PrePopulatePool(Pool* pool, GsampleInfo* sampleInfo);
void MutateNeRnd(Bicho* bicho, Pool* pool);
void MutateLateral(Bicho* bicho, Pool* pool);
void Run(Pool* pool, GsampleInfo* sampleInfo, std::string fichSal);
void RunDbg(Pool* pool, GsampleInfo* sampleInfo, std::string fichSal);

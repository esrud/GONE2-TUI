#pragma once

// Runtime knobs for the GA loop. Each `Pool` carries one of these
// and the smoothness branch in `CalculaSCImpl` reads the matching
// copy from `GsampleInfo`. The values default to whatever the
// compile-time -DGA_* flags set at build time, so non-combo builds
// behave bit-identically to the pre-runtime-config code path; the
// `gone-ncurses_combo` target swaps the struct between three named
// presets without re-reading input or recomputing d².

enum class GASmoothKind {
  Quadratic = 0,  // smooth += d * d           (L2)
  L1        = 1,  // smooth += |d|             (total variation)
  Truncated = 2,  // smooth += min(d², cutoff²) (truncated quadratic)
};

struct GAConfig {
  double smoothLambda;
  GASmoothKind smoothKind;
  double smoothCutoff;
  double kickProb;
  double kickMag;
  double scoutFrac;
  // Heterozygosity anchor on the deep-time end of the Ne curve.
  // CalculaSC adds  λ_anchor · (log10(Necons) − log10(anchorNeHaploid))²
  // to the score so the GA's ancient plateau is softly pulled toward
  // the mutation-drift-balance target H/(4μ). Both fields default
  // to 0 (anchor off); libgone fills them when -m <μ> is given.
  // anchorNeHaploid is in haploid units (matches Bicho::NeBl).
  double anchorNeHaploid;
  double anchorLambda;
  // Short label shown in the live TUI status while this config runs
  // (e.g. "trunc05_kick"). Empty in non-combo builds.
  const char* label;
};

// Build the default config from the compile-time -DGA_* flags that
// the rest of the codebase was built with. Used both as the
// configuration in non-combo builds and as one of the three combo
// presets.
inline GAConfig MakeDefaultGAConfig() {
  GAConfig c{};
#ifdef GA_SMOOTH_LAMBDA
  c.smoothLambda = static_cast<double>(GA_SMOOTH_LAMBDA);
#else
  c.smoothLambda = 0.0001;
#endif
#if defined(GA_SMOOTH_L1)
  c.smoothKind = GASmoothKind::L1;
#elif defined(GA_SMOOTH_CUTOFF)
  c.smoothKind = GASmoothKind::Truncated;
#else
  c.smoothKind = GASmoothKind::Quadratic;
#endif
#ifdef GA_SMOOTH_CUTOFF
  c.smoothCutoff = static_cast<double>(GA_SMOOTH_CUTOFF);
#else
  c.smoothCutoff = 0.5;
#endif
#ifdef GA_KICK_PROB
  c.kickProb = static_cast<double>(GA_KICK_PROB);
#else
  c.kickProb = 0.05;
#endif
#ifdef GA_KICK_MAG
  c.kickMag = static_cast<double>(GA_KICK_MAG);
#else
  c.kickMag = 10.0;
#endif
#ifdef GA_SCOUT_FRAC
  c.scoutFrac = static_cast<double>(GA_SCOUT_FRAC);
#else
  c.scoutFrac = 0.4;
#endif
  c.anchorNeHaploid = 0.0;
  c.anchorLambda    = 0.0;
  c.label = "";
  return c;
}

// The three presets used by the combo target. Order is the order
// of the sweep loop in libgone.cpp.
inline GAConfig MakeComboTruncKickConfig() {
  GAConfig c = MakeDefaultGAConfig();
  c.smoothLambda = 0.0001;
  c.smoothKind   = GASmoothKind::Truncated;
  c.smoothCutoff = 0.5;
  c.kickProb     = 0.15;
  c.kickMag      = 30.0;
  c.scoutFrac    = 0.7;
  c.label        = "trunc05_kick";
  return c;
}

inline GAConfig MakeComboL2Config() {
  GAConfig c = MakeDefaultGAConfig();
  c.smoothLambda = 0.0001;
  c.smoothKind   = GASmoothKind::Quadratic;
  c.smoothCutoff = 0.5;
  c.kickProb     = 0.05;
  c.kickMag      = 10.0;
  c.scoutFrac    = 0.4;
  c.label        = "l2";
  return c;
}

inline GAConfig MakeComboL1KickConfig() {
  GAConfig c = MakeDefaultGAConfig();
  c.smoothLambda = 0.0001;
  c.smoothKind   = GASmoothKind::L1;
  c.smoothCutoff = 0.5;
  c.kickProb     = 0.15;
  c.kickMag      = 30.0;
  c.scoutFrac    = 0.7;
  c.label        = "l1_kick";
  return c;
}

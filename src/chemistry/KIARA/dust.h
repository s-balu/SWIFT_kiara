/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Balu Sreedhar
 *
 * Compile-time layout and small helpers for the KIARA binned dust model.
 ******************************************************************************/
#ifndef SWIFT_KIARA_DUST_H
#define SWIFT_KIARA_DUST_H

/* config.h normally provides this.  Keep headers usable by stand-alone tools. */
#ifndef KIARA_DUST_N_BINS
#define KIARA_DUST_N_BINS 6
#endif

#define KIARA_DUST_N_MATERIALS 2
#define KIARA_DUST_N_DISTRIBUTION_VALUES \
  (KIARA_DUST_N_MATERIALS * KIARA_DUST_N_BINS)

enum kiara_dust_material {
  kiara_dust_carbonaceous = 0,
  kiara_dust_silicate = 1,
  kiara_dust_material_count = KIARA_DUST_N_MATERIALS
};

/** Zero a material-by-size distribution. */
static inline void
kiara_dust_distribution_zero(
    float distribution[KIARA_DUST_N_MATERIALS][KIARA_DUST_N_BINS]) {
  for (int material = 0; material < KIARA_DUST_N_MATERIALS; ++material)
    for (int bin = 0; bin < KIARA_DUST_N_BINS; ++bin)
      distribution[material][bin] = 0.f;
}

/** Return the sum of all material-by-size mass fractions. */
static inline float
kiara_dust_distribution_sum(
    const float distribution[KIARA_DUST_N_MATERIALS][KIARA_DUST_N_BINS]) {
  float sum = 0.f;
  for (int material = 0; material < KIARA_DUST_N_MATERIALS; ++material)
    for (int bin = 0; bin < KIARA_DUST_N_BINS; ++bin)
      sum += distribution[material][bin];
  return sum;
}

/** Clamp negative entries and normalize a non-empty distribution to unity. */
static inline void
kiara_dust_distribution_normalize(
    float distribution[KIARA_DUST_N_MATERIALS][KIARA_DUST_N_BINS]) {
  float sum = 0.f;
  for (int material = 0; material < KIARA_DUST_N_MATERIALS; ++material) {
    for (int bin = 0; bin < KIARA_DUST_N_BINS; ++bin) {
      if (distribution[material][bin] < 0.f)
        distribution[material][bin] = 0.f;
      sum += distribution[material][bin];
    }
  }
  if (sum > 0.f) {
    const float sum_inv = 1.f / sum;
    for (int material = 0; material < KIARA_DUST_N_MATERIALS; ++material)
      for (int bin = 0; bin < KIARA_DUST_N_BINS; ++bin)
        distribution[material][bin] *= sum_inv;
  }
}

/* Assign a smooth default mass distribution over logarithmically spaced bins.
 * The caller may overwrite this with a model-specific source distribution. */
static inline void kiara_dust_distribution_powerlaw(
    float distribution[KIARA_DUST_N_MATERIALS][KIARA_DUST_N_BINS]) {
  kiara_dust_distribution_zero(distribution);
  for (int material = 0; material < KIARA_DUST_N_MATERIALS; ++material) {
    for (int bin = 0; bin < KIARA_DUST_N_BINS; ++bin)
      distribution[material][bin] = 1.f / (float)KIARA_DUST_N_BINS;
  }
  kiara_dust_distribution_normalize(distribution);
}

#endif /* SWIFT_KIARA_DUST_H */

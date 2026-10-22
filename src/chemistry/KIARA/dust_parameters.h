/* Shared runtime grain-radius parameters for KIARA cooling and feedback. */
#ifndef SWIFT_KIARA_DUST_PARAMETERS_H
#define SWIFT_KIARA_DUST_PARAMETERS_H

#include "dust_grid.h"
#include "error.h"
#include "parser.h"

/** Read one common grid, retaining the legacy endpoint parameter aliases. */
static inline void kiara_dust_read_grid(struct swift_params *params,
                                        double radius[KIARA_DUST_N_BINS]) {
  const double legacy_min = parser_get_opt_param_double(
      params, "KIARACooling:dust_small_grainsize", 0.005);
  const double legacy_max = parser_get_opt_param_double(
      params, "KIARACooling:dust_large_grainsize", 0.1);
  const double minimum = parser_get_opt_param_double(
      params, "KIARACooling:dust_grain_size_min_micron", legacy_min);
  const double maximum = parser_get_opt_param_double(
      params, "KIARACooling:dust_grain_size_max_micron", legacy_max);
  if (parser_does_param_exist(params, "KIARACooling:dust_grain_sizes_micron")) {
    parser_get_opt_param_double_array(params,
                                      "KIARACooling:dust_grain_sizes_micron",
                                      KIARA_DUST_N_BINS, radius);
  } else if (!kiara_dust_grid_logarithmic(minimum, maximum, radius)) {
    error("KIARA dust grid requires finite radii with 0 < minimum < maximum");
  }
  if (!kiara_dust_grid_valid(radius))
    error(
        "KIARACooling:dust_grain_sizes_micron must be finite, positive and "
        "strictly increasing");
}

#endif

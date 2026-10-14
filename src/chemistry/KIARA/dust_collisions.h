/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Balu Sreedhar
 *
 * Shared numerical helpers for KIARA grain collisions.
 ******************************************************************************/
#ifndef SWIFT_KIARA_DUST_COLLISIONS_H
#define SWIFT_KIARA_DUST_COLLISIONS_H

#include "dust.h"

#include <assert.h>
#include <math.h>

/** Representative radii and single-grain masses, in cgs units.
 *
 * The existing KIARA radii are mass pivots, not bin edges. Collision products
 * are interpolated between these pivots to avoid suppressing coagulation on
 * the coarse six-bin grid. No extra particle fields or size bins are needed.
 */
struct kiara_dust_collision_grid {
  double radius[KIARA_DUST_N_BINS];
  double mass[KIARA_DUST_N_BINS];
};

/**
 * @brief Build a collision grid from radii in micron and density in g/cm^3.
 * @return One on success, zero for non-finite, non-positive or unordered data.
 * The output is only changed on success.
 */
static inline int kiara_dust_collision_grid_init(
    struct kiara_dust_collision_grid *grid,
    const double radius_micron[KIARA_DUST_N_BINS], const double grain_density) {

  if (!isfinite(grain_density) || grain_density <= 0.) return 0;

  struct kiara_dust_collision_grid candidate;
  const double four_pi_over_three = 4.1887902047863909846;
  for (int i = 0; i < KIARA_DUST_N_BINS; ++i) {
    candidate.radius[i] = radius_micron[i] * 1.e-4;
    const double a = candidate.radius[i];
    candidate.mass[i] = four_pi_over_three * grain_density * a * a * a;
    if (!isfinite(a) || a <= 0. || !isfinite(candidate.mass[i]) ||
        candidate.mass[i] <= 0.)
      return 0;
    if (i > 0 && (a <= candidate.radius[i - 1] ||
                  candidate.mass[i] <= candidate.mass[i - 1]))
      return 0;
  }
  *grid = candidate;
  return 1;
}

/**
 * @brief Deposit collision products into one material's mass array.
 *
 * @param grid Valid collision grid.
 * @param grain_mass Mass of one product grain, in grams (positive and finite).
 * @param total_mass Total mass of all such products (non-negative and finite).
 * May also be a mass fraction: the destination uses the same units.
 * @param bin_mass Destination material row, updated additively.
 *
 * Between pivots, interpolate product NUMBER linearly in grain mass, then
 * convert to bin masses. This conserves both deposited mass and grain number:
 * sum(delta_M_i) = total_mass and sum(delta_M_i / m_i) = total_mass/grain_mass.
 * Interpolating mass fractions directly would violate the number constraint.
 *
 * Below/above the represented range, retain all mass in the nearest endpoint.
 * This closed-boundary approximation conserves material mass, but cannot
 * conserve grain number outside the grid. It must not be interpreted as
 * vaporization or as resolved growth beyond the largest radius.
 *
 * This is a coarse-grid adaptation of Li et al. (2021), section 2.5.3;
 * their 41-bin, piecewise-constant discretization is not used here.
 */
static inline void kiara_dust_collision_deposit(
    const struct kiara_dust_collision_grid *grid, const double grain_mass,
    const double total_mass, double bin_mass[KIARA_DUST_N_BINS]) {

  assert(isfinite(grain_mass) && grain_mass > 0.);
  assert(isfinite(total_mass) && total_mass >= 0.);
  if (total_mass == 0.) return;

  if (grain_mass <= grid->mass[0]) {
    bin_mass[0] += total_mass;
    return;
  }
  const int last = KIARA_DUST_N_BINS - 1;
  if (grain_mass >= grid->mass[last]) {
    bin_mass[last] += total_mass;
    return;
  }

  for (int hi = 1; hi < KIARA_DUST_N_BINS; ++hi) {
    if (grain_mass <= grid->mass[hi]) {
      const double lo_mass = grid->mass[hi - 1];
      const double hi_mass = grid->mass[hi];
      const double lo_fraction = (lo_mass / grain_mass) *
                                 ((hi_mass - grain_mass) / (hi_mass - lo_mass));
      const double deposited_lo = total_mass * lo_fraction;
      bin_mass[hi - 1] += deposited_lo;
      bin_mass[hi] += total_mass - deposited_lo;
      return;
    }
  }
}

/**
 * @brief Li et al. (2021), equation 15, grain velocity dispersion in cm/s.
 *
 * All inputs are cgs: grain radius in cm, gas/grain densities in g/cm^3,
 * sound speed in cm/s. This returns the THREE-dimensional dispersion sigma_gr;
 * each Gaussian velocity component has standard deviation sigma_gr/sqrt(3).
 * The caller chooses the physical gas phase and performs code-unit conversion.
 *
 * L_J = sqrt(pi c_g^2 / (G rho_g))/2 and v_J = 0.7 sqrt(L_J/pc) km/s.
 * Reference: https://arxiv.org/abs/2012.03978, section 2.5.3.
 */
static inline double kiara_dust_collision_velocity_dispersion(
    const double radius, const double gas_density, const double sound_speed,
    const double grain_density) {

  assert(isfinite(radius) && radius > 0.);
  assert(isfinite(gas_density) && gas_density > 0.);
  assert(isfinite(sound_speed) && sound_speed > 0.);
  assert(isfinite(grain_density) && grain_density > 0.);

  const double pi = 3.14159265358979323846;
  const double gravity_cgs = 6.67430e-8;
  const double hydrogen_mass_cgs = 1.6735575e-24;
  const double parsec_cgs = 3.0856775814913673e18;
  const double jeans_length =
      0.5 * sound_speed * sqrt(pi / (gravity_cgs * gas_density));
  const double mach_squared =
      (7.e4 / sound_speed) * (7.e4 / sound_speed) * (jeans_length / parsec_cgs);

  return 6.e3 * mach_squared * (radius / 1.e-5) *
         sqrt(hydrogen_mass_cgs / gas_density) * (grain_density / 2.4);
}

#endif /* SWIFT_KIARA_DUST_COLLISIONS_H */

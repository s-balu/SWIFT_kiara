KIARA dust
==========

Grain-size grid
--------------

The gas-carried dust model tracks carbonaceous and silicate mass on a common
radius grid. Choose the number of radii per material at build time with
``--with-kiara-dust-size-bins=N`` (default 6, minimum 2). The number sets the
particle layout; changing it requires rebuilding and does not convert existing
restart files.

Set ``KIARACooling:dust_grain_size_min_micron`` and
``KIARACooling:dust_grain_size_max_micron`` at run time to construct logarithmic
representative radii including both endpoints. The defaults remain 0.005 and
0.1 micron. For Li et al. (2021)'s wider range, set these to 0.0001 and 1.
The legacy keys ``dust_small_grainsize`` and ``dust_large_grainsize`` remain
fallback aliases when the corresponding new key is absent.

Alternatively, ``KIARACooling:dust_grain_sizes_micron`` supplies exactly N
finite, positive, strictly increasing radii and overrides the endpoint keys.
An explicit grid may be nonuniform. The radii are physical grain radii, not
comoving lengths, and are shared by cooling and stellar dust injection.

Surviving dust below the minimum represented radius is deposited in the
smallest bin; dust above the maximum is deposited in the largest bin. This
retains mass at the boundaries but cannot retain the unresolved grain number
or surface area. Physically destroyed dust is returned to gas-phase elements.
Six bins across a wider range remain a coarse approximation; configurable
geometry does not imply convergence with respect to resolution.

``KIARACooling:dust_grainsize`` is the separate representative radius passed
to Grackle and does not define the framework's size grid.

Dust sources
------------

``KIARACooling:dust_initial_source`` selects the Li source distribution for
dust without an existing size distribution, empty material rows, and the
artificial self-enrichment seed. It defaults to ``SNII`` and also accepts
``AGB``. This is an explicit assumption about dust without a recorded source
history, not an inference of that history from its composition. The source
parameters and integration are shared with stellar injection; there is no
separate single-radius initialization prescription.

Stellar injection integrates the mass spectrum corresponding to Li et al.
(2021), equation 5, over geometric bin cells. The decaying lognormal form is
used (the printed equation omits the minus sign in the exponent). Both tails
are added to the endpoint bins rather than discarded or sampled randomly.
The number spectrum is proportional to
``a**(-p) * exp(-ln(a/a0)**2/(2*sigma**2))`` per logarithmic radius interval;
mass weighting contributes another factor of ``a**3``.

The source parameters are ``KIARAFeedback:dust_SNII_a0_micron`` (0.1),
``dust_SNII_sigma`` (0.6), and ``dust_SNII_p`` (0), with corresponding
``dust_AGB_a0_micron`` (0.1), ``dust_AGB_sigma`` (0.47), and ``dust_AGB_p`` (4).
The old ``dust_small_fraction_SNII/AGB`` feedback inputs are obsolete and
produce a warning if supplied. Both materials use the same source spectrum
for a given channel; their masses and subsequent size distributions remain
independent. The existing 100 Myr stellar-age channel classification remains
an approximation; channel separation using SSP tables is unchanged.

Radius evolution
----------------

Growth applies ``a_new = a + a_ref * dt / tau_base`` from Li et al. (2021),
equations 6--7, with ``a_ref = 0.1 micron`` a physical reference rather than
a grid endpoint. Each grain's mass changes as ``a_new**3`` before conservative
deposition. Only the size shape within each material is updated here: Grackle
continues to own the accreted mass and elemental depletion, using its separate
representative grain radius. This remains an approximate coupling rather than
a fully size-resolved gas-metal accretion calculation.

Thermal sputtering applies the analytic radius decrement implied by equations
11--12 under the local conditions held fixed for the step. Surviving mass is
computed from the remaining grain volume in each bin. Carbon and non-carbon
elemental losses are returned to the corresponding gas inventories; the
previous total-mass destruction cap and adjacent-bin shuffle are removed.
Positive sub-grid radii retain their surviving mass in the smallest bin;
zero-radius grains have been physically destroyed. Each transport pass uses
the old distribution once, so grains cannot be processed repeatedly while
looping over newly populated bins.

Validation
----------

Run ``python3 tests/testKIARADustGrid.py`` from the repository root to compile
standalone checks with 2, 3, 6, 12 and 41 bins under address and undefined
behaviour sanitizers. The checks include nonuniform grids, analytic source
integrals, physical-radius remapping, collision boundary retention, and the
actual cooling dust routines with controlled hydro/unit adapters at scale
factors 1, 0.5 and 0.1. These do not link Grackle or GSL and do not replace
a full configured KIARA build or an evolution convergence study.

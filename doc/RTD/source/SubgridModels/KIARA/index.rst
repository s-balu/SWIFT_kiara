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

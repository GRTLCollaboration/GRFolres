# KerrBH4dST

A single Kerr black hole in 4-derivative scalar-tensor gravity, ported from
`GRChombo/Examples/KerrBH4dST/` of GRFolres to GRTeclyn.

Evolved variables: the CCZ4 variables plus the scalar field `phi` and its
conjugate momentum `Pi` (`StateVariables.hpp`).

## Initial data

- Metric: GRTeclyn's `KerrBHInitialData` (`Source/BlackHoles/`), the
  quasi-isotropic Kerr solution of arXiv:1001.4077, with a pre-collapsed
  lapse `sqrt(chi)`. Parameters: `kerr.mass`, `kerr.spin` (a = J/M,
  |a| <= M), optional `kerr.spin_direction` and `kerr.offset`.
- Scalar field: a Gaussian `phi` profile (`ScalarFieldInitialData.hpp`,
  `scalar_field.*` parameters) with `Pi = 0`.
- `Gamma^i` is computed numerically from the metric (`GammaCalculator`), and
  `B^i` is set to `shift_Gamma_coeff * Gamma^i - eta * shift^i`, so that the
  nonzero Kerr shift is initially approximately stationary under the
  modified puncture gauge (`B^i` is not evolved). This differs from the
  GRChombo version, where `B^i = 0`.

## Symmetry

A spinning BH is only symmetric under z -> -z, so the parameter files use a
half domain with a reflective lower z boundary. The spin must therefore be
along z. Use Sommerfeld boundaries everywhere (full domain) for any other
spin direction.

## Building and running

```bash
make -j GRTECLYN_HOME=/path/to/GRTeclyn AMREX_HOME=/path/to/amrex
./KerrBH4dST3d.<comp>.ex params.txt        # example run
./KerrBH4dST3d.<comp>.ex params_test.txt   # 8-step quick test (~10 s on CPU)
```

`params.txt` uses a coarse two-level fixed grid that does not resolve the
horizon well; increase `amr.max_level` for production runs.

`params.txt` writes `phi` along the diagonal to
`extraction_data/phi_profile_*.dat`. `plot_lineouts.py` plots `phi` and `rho`
lineouts from `data/`; it needs the `rho` extraction, which is not yet
implemented, so adjust it (or `data_path`) before use.

## Not yet ported

- Constraints, Weyl4 and rho diagnostics (commented out in
  `KerrBH4dSTLevel.cpp`).
- Apparent-horizon finder.

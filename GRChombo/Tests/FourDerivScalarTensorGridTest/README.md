# FourDerivScalarTensorGridTest (GRChombo)

Reference-value generator for the GRTeclyn test
`GRFolres/GRTeclyn/Tests/FourDerivScalarTensorTest`.

Unlike `../FourDerivScalarTensorTest` (which feeds hand-built derivatives into
the RHS functions and checks against Mathematica), this program fills a grid
with the same polynomial initial data as the GRTeclyn test, runs
`ModifiedCCZ4RHS::compute` (vacuum CCZ4 + modified gauge + effective EM tensor +
scalar evolution + principal-part solve + KO dissipation), and prints the RHS at
the probe cell as a block of

```
known[c_chi] = ...;
...
have_reference = true;
```

ready to paste into
`GRFolres/GRTeclyn/Tests/FourDerivScalarTensorTest/values1.hpp`.

## Build & run

Needs `CHOMBO_HOME` set and a Chombo build (same as the other GRChombo tests):

```bash
cd GRFolres/GRChombo/Tests/FourDerivScalarTensorGridTest
make all
./FourDerivScalarTensorGridTest*.ex
```

## Keeping the two sides in sync

The grid size, `dx`, probe cell, CCZ4 / gauge / modified-gauge parameters and the
coupling constants must match
`GRFolres/GRTeclyn/Tests/FourDerivScalarTensorTest/` exactly:

| quantity            | value                              |
|---------------------|------------------------------------|
| `N_GRID`            | 16                                 |
| ghost cells         | 3                                  |
| `dx`                | `0.5 / N_GRID`                      |
| probe cell          | `(8, 8, 8)`                         |
| `x` at index `i`    | `i * dx`  (node convention)        |
| `kappa1,2,3`        | `0.1, 0.0, 1.0`                     |
| `covariantZ4`       | `0`                                |
| `formulation`       | CCZ4                               |
| `sigma`             | `1.0`                              |
| lapse coeffs        | advec `1`, power `1`, coeff `2`     |
| shift coeffs        | advec `1`, Gamma `0.75`, eta `1`    |
| `a0`, `b0`          | `0.35`, `0.55`                      |
| `G_Newton`          | `1.0`                              |
| coupling constants  | see `CouplingAndPotential.hpp`      |

If the 4dST equations or these constants change, rebuild, rerun, and repaste
`values1.hpp`.

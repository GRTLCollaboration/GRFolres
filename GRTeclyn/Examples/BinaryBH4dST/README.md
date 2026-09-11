# BinaryBH4dST

Binary black holes in shift-symmetric 4-derivative scalar-tensor gravity
(Einstein-scalar-Gauss-Bonnet), ported from
`GRChombo/Examples/BinaryBH4dST/` of GRFolres to GRTeclyn.

Evolved variables: the CCZ4 variables plus the scalar field `phi` and its
conjugate momentum `Pi` (`StateVariables.hpp`).

Diagnostics (AMReX derived records, computed on demand for plotfiles /
extraction):

| record             | components                             | class |
|--------------------|----------------------------------------|-------|
| `constraints`      | `Ham`, `Mom1`, `Mom2`, `Mom3`          | `ModifiedGravityConstraints` |
| `Weyl4`            | `Weyl4_Re`, `Weyl4_Im`                | `ModifiedGravityWeyl4` |
| `rho_diagnostics`  | `rho_phi`, `rho_g2`, `rho_g3`, `rho_GB` | `RhoDiagnostics` |

## Building

The example depends on the modified-gravity source that lives in
`GRTeclyn/Source/{FourDerivScalarTensor,ModifiedCCZ4}/` (ported separately).
Assumed interface points are listed in the header comment of
`BinaryBH4dSTLevel.cpp`; reconcile the names/signatures there with the Source
port when it lands.

```bash
make -j                       # boosted-BH superposed initial data
make -j USE_TWOPUNCTURES=TRUE TWOPUNCTURES_SOURCE=/path/to/TwoPunctures/Source
```

## Not yet ported

- Constraint-norm output (`calculate_constraint_norms`) — needs an
  `AMRReductions` equivalent in GRTeclyn.
- Apparent-horizon finder (`AH_*` parameters).

Use `KerrBH4dST/` for a fixed-grid single-BH 4dST example.

claude --resume 9b892624-3162-47ef-a626-c02444f69a42


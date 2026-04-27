# Results Directory

## Overview
This folder contains outputs from the WCVI-E2E model simulations, including long-term simulations generated during PhD research.

---

## Contents

### `outputs_physics/` or `outputs_e2e/`
These folders are created when the physics-only version and the full E2E version are run. They will contain the archived model outputs. Refer to scripts/WCVI-E2E_run.m

---

### `100-year simulations/`
Long-term (e.g., 100-year) E2E simulations used in the PhD analysis, with each simulation corresponding to a calibration step (see docs/thesis or docs/paper_submitted)
1. Baselinesim = baseline simulation 
2. Linear = adjusted trophic parameters and linear non-predatory mortality formulation
3. QuadraticPLANKTON = adjusted trophic parameters and quadratic non-predatory mortality formulation for plankton groups (and linear non-predatory mortality formulation for nektonic groups)
4. Quadratic = adjusted trophic parameters and quadratic non-predatory mortality formulation for all groups

**Note:**
- These files may be computationally expensive to reproduce.
- They are included for reference and comparison with published results.
- Only the final year of 100-year simulations is included, where the system was considered to have reached equilibrium. Full 100-year simulation outputs are excluded due to file size constraints.

---

## File Types
- Tables
- Monthly time-series data of ecosystem state variables, and diagnostic variables (derivatives, intermediate fluxes, additional diagnostic variables)   

Refer to documentation in the `scripts/` directory for variable definitions and analysis workflows.

---

## Reproducing Results
To generate new simulations:

1. physics-only version

```
physicalmodel('outputs_physics');
```

2. Full E2E version

```
physicalmodel('outputs_e2e','biofun',@biomodel);
```

3. Outputs will be saved to `outputs_physics/` or `outputs_e2e/`

**Note:**
- Results are specific to the WCVI ecosystem configuration used in this study

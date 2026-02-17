# unresolvedSphereVoFPFPlus / fallingSphere

This tutorial is an **unresolved VOF-DEM** case for testing the new
`unresolvedSphereVoFPFPlus` solver.

## What it does
- Uses a two-phase VOF setup (`alpha.water` + `p_rgh`) similar to interFoam.
- Uses unresolved DEM coupling (`constant/couplingProperties` -> `unresolved`).
- Runs CFD in parallel with `unresolvedSphereVoFPFPlus`.

## Run
From this directory:

```bash
./Allrun
```

## Expected main outputs
- OpenFOAM time folders with `U`, `p_rgh`, `alpha.water`, etc.
- DEM coupling data written by phasicFlow.
- VTK outputs from:
  - `foamToVTK -time 0:0.5`
  - `pFlowToVTK -t 0:0.5`

## Notes
- Requires OpenFOAM environment + phasicFlowPlus tools in PATH.
- Solver executable must be compiled first (`unresolvedSphereVoFPFPlus`).

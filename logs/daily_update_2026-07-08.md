# Daily update - 2026-07-08

## Heat-conduction case

- Added `examples/complex_heat_multimaterial.plc`, a larger 13-vertex concave footprint with extrusion height `1.60`.
- Added five material regions over a default matrix material:
  - copper core (`k = 385.0`)
  - ceramic cap (`k = 24.0`)
  - insulation band (`k = 0.08`)
  - graphite midplane (`k = 130.0`)
  - aluminum lip (`k = 205.0`)

## Code improvements

- Extended PLC parsing to skip comments and read optional default conductivity plus material-region records.
- Added per-wedge conductivity/material metadata through extrusion, OpenVolumeMesh round-trip, partition reordering, and FEM assembly.
- Improved arbitrary polygon filtering by keeping triangle centroids inside the outer polygon instead of relying on a single seed triangle.
- Added mesh boundary metadata so Dirichlet vertices are detected correctly for non-rectangular footprints.

## Verification

- Pending build and MPI run after the implementation commit, per branch workflow.

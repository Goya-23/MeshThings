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

- Configured and built successfully:
  - `CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release`
  - `cmake --build build -j`
- Ran the new complex multi-material case:
  - `mpirun --allow-run-as-root -np 2 ./build/heat_conduction_mpi examples/complex_heat_multimaterial.plc i 18 18 6`
  - Mesh: 302 footprint vertices, 451 triangles, 2114 wedge vertices, 2706 wedge cells.
  - Material counts: default 1147, copper core 378, ceramic cap 458, insulation band 258, graphite midplane 312, aluminum lip 153 wedge cells.
  - HYPRE PCG converged in 24 iterations with residual `2.7354e-11`.
- Re-ran the original unit brick path for backward compatibility:
  - `mpirun --allow-run-as-root -np 2 ./build/heat_conduction_mpi examples/unit_brick.plc i 6 6 3`
  - HYPRE PCG converged in 6 iterations with residual `1.72381e-12`.

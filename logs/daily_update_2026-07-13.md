# Daily update - 2026-07-13

## Large heat-conduction benchmark case

- Added `examples/large_composite_heat_multimaterial.plc`, a 20-vertex concave composite footprint extruded to height `2.40`.
- Added twelve material regions over a default matrix material, including copper, graphite, ceramic, polymer, steel, aluminum, titanium, foam, silicon, air, silver, and grease regions.
- Target verification command:
  - `mpirun --allow-run-as-root -np <N> ./build/heat_conduction_mpi examples/large_composite_heat_multimaterial.plc i 110 90 10`
- Verified mesh density exceeds the requested threshold with 106,230 wedge cells.

## Code quality and benchmark instrumentation

- Reduced FEM assembly work by assembling only wedges incident to vertices owned by the current MPI rank.
- Added benchmark telemetry to `examples/heat_conduction_mpi.cpp`:
  - per-stage min/average/max timings across ranks
  - local sparse matrix row and nonzero balance
  - RSS and peak RSS memory summaries
  - machine-parseable `BENCHMARK` output

## Verification and process-count comparison

- Installed missing runtime/build packages in this environment:
  - `openmpi-bin libopenmpi-dev libmetis-dev libhypre-dev time`
- Configured and built successfully:
  - `CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release`
  - `cmake --build build -j`
- Official benchmark mesh:
  - Footprint CDT: 9,721 vertices, 10,623 triangles
  - Wedge mesh: 106,931 vertices, 106,230 wedge cells
  - Matrix: 106,931 rows, 992,274 nonzeros
  - HYPRE PCG iterations: 539
  - Relative L2 error: `0.355407`

| MPI processes | Total seconds | Peak RSS max (MB) | Assembly seconds max | HYPRE solve seconds max | Raw log |
| ---: | ---: | ---: | ---: | ---: | --- |
| 1 | 4.95674 | 277.492 | 0.511282 | 3.35609 | `logs/benchmarks/large_composite_np1_2026-07-13.log` |
| 2 | 4.80963 | 222.293 | 0.277842 | 3.38183 | `logs/benchmarks/large_composite_np2_2026-07-13.log` |
| 4 | 4.74083 | 218.25 | 0.183385 | 3.41238 | `logs/benchmarks/large_composite_np4_2026-07-13.log` |

Notes:

- `110 90 8` was tried first and produced 84,984 wedge cells, so the official run was increased to `nz = 10`.
- Assembly time improved with more ranks because each rank now assembles only its incident wedges.
- Total runtime remains dominated by the current gather-to-rank-0 HYPRE solve path.

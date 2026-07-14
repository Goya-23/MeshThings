# Daily update - 2026-07-14

## Heat-conduction case

- Added `examples/large_composite_heat_multimaterial.plc`, a 22-vertex
  non-rectangular composite footprint extruded to height `2.40`.
- Added 12 material regions over a default matrix material, covering high
  conductivity copper/diamond spreaders, moderate metal/ceramic/graphite
  regions, and low conductivity insulation/polymer/coolant regions.
- Target benchmark resolution is `nx = 110`, `ny = 90`, `nz = 10`, which is
  intended to exercise more than 100,000 wedge cells.

## Code quality and instrumentation

- Updated MPI heat-driver output with stage timings, row/nnz/incident-wedge
  balance, peak RSS memory summaries, and a machine-parseable `BENCHMARK` line.
- Scoped FEM assembly to wedges incident to each rank's owned vertices instead
  of assembling the full global sparse map on every rank.

## Verification

- Installed required native dependencies in the cloud environment:
  `openmpi-bin`, `libopenmpi-dev`, `libmetis-dev`, `libhypre-dev`, and `time`.
- Configured and built successfully:
  - `CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release`
  - `cmake --build build -j`
- Ran the unit-brick MPI regression:
  - `mpirun --allow-run-as-root -np 2 ./build/heat_conduction_mpi examples/unit_brick.plc i 6 6 3`
  - HYPRE PCG converged in 6 iterations with residual `1.72381e-12`.
- Ran the new large composite case at `nx = 110`, `ny = 90`, `nz = 10`:
  - Mesh: 9,723 footprint vertices, 11,065 triangles, 106,953 wedge vertices,
    110,650 wedge cells.
  - HYPRE PCG converged in 427-428 iterations depending on process count.
  - Relative L2 error vs manufactured solution: `0.385714`.

## Process-count benchmark

| MPI processes | Driver total (s) | `/usr/bin/time` wall (s) | Assembly (s) | Solve (s) | Peak RSS max (KB) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 4.7185 | 5.12 | 0.616573 | 2.8574 | 291728 |
| 2 | 4.5285 | 4.92 | 0.369824 | 2.83551 | 232888 |
| 4 | 4.92442 | 5.35 | 0.294131 | 3.25617 | 229872 |

- Assembly time improved with more MPI ranks after scoping assembly to incident
  wedges.
- End-to-end runtime did not scale monotonically because mesh generation still
  happens independently on every rank and HYPRE still solves the gathered matrix
  serially on rank 0.
- Detailed benchmark records are in
  `logs/benchmarks/large_composite_heat_multimaterial_2026-07-14.md` and the
  corresponding `_np1.txt`, `_np2.txt`, and `_np4.txt` raw logs.

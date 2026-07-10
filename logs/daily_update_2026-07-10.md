# Daily update - 2026-07-10

## Large heat-conduction case

- Added `examples/large_composite_heat_multimaterial.plc`.
- The case uses a 24-vertex concave 2D footprint with extrusion height `2.40`.
- Added twelve material regions over the default composite matrix material (`k = 2.5`):
  - copper bus (`k = 385.0`)
  - aluminum frame (`k = 205.0`)
  - ceramic plate (`k = 24.0`)
  - graphite spreader (`k = 130.0`)
  - insulation slot (`k = 0.06`)
  - steel anchor (`k = 45.0`)
  - silicon die (`k = 148.0`)
  - polymer lid (`k = 0.22`)
  - coolant channel (`k = 0.60`)
  - thermal interface (`k = 8.0`)
  - tungsten insert (`k = 170.0`)
  - air gap (`k = 0.026`)

## Code improvements

- Localized heat FEM assembly so each MPI rank only assembles wedges incident to its owned vertices instead of assembling the full mesh on every rank.
- Added stage timing, sparse matrix row/nnz balance, and RSS/peak-RSS reporting to `examples/heat_conduction_mpi.cpp`.
- Added a machine-parseable `BENCHMARK` line to each run for reproducible performance logs.

## Verification and benchmarks

- Configured and built successfully:
  - `CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release`
  - `cmake --build build -j`
- Smoke-tested the original unit brick path:
  - `mpirun --allow-run-as-root -np 2 ./build/heat_conduction_mpi examples/unit_brick.plc i 6 6 3`
- Validated the new composite geometry/material setup:
  - `mpirun --allow-run-as-root -np 2 ./build/heat_conduction_mpi examples/large_composite_heat_multimaterial.plc i 30 25 3`
- Ran the large `110 90 9` mesh with 1, 2, and 4 MPI ranks:
  - Mesh: 9,725 footprint vertices, 11,725 footprint triangles, 97,250 wedge vertices, 105,525 wedge cells.
  - Matrix: 942,974 nonzeros.
  - HYPRE PCG: 526 iterations, relative L2 error `0.395812`.

| MPI ranks | Internal total (s) | Assembly (s) | Solve (s) | Peak RSS/rank (MB) | External wall time |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 4.757378 | 0.519714 | 3.135757 | 268.804688 | 0:05.13 |
| 2 | 4.570183 | 0.285066 | 3.122220 | 215.320312 | 0:04.93 |
| 4 | 4.568784 | 0.183612 | 3.158071 | 210.800781 | 0:04.95 |

Detailed benchmark logs are in `logs/benchmarks/`.

## Observations

- The new large case satisfies the `>100000` mesh-element target with 105,525 wedge cells.
- Assembly runtime and per-rank memory improve with more MPI ranks after the assembly locality change.
- Overall runtime is dominated by the current gather-to-rank-0 solver path, so adding ranks primarily improves assembly and per-rank memory rather than total solve time.

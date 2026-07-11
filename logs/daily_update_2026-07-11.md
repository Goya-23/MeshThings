# Daily update - 2026-07-11

## Large heat-conduction case

- Added `examples/large_composite_heat_multimaterial.plc`.
- The case uses a 24-vertex concave 2D footprint with extrusion height `2.40`.
- The benchmarked `110 90 9` discretization creates 105,525 wedge cells, satisfying the `>100000` mesh-element target.
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

## Code improvements and debugging

- Localized heat FEM assembly so each MPI rank only assembles wedges incident to its owned vertices instead of assembling the full mesh on every rank.
- Added stage timing, sparse matrix row/nnz balance, and RSS/peak-RSS reporting to `examples/heat_conduction_mpi.cpp`.
- Added a machine-parseable `BENCHMARK` line to each run for reproducible performance logs.
- Initial configure failed in this environment because MPI was missing; installed `openmpi-bin`, `libopenmpi-dev`, `libmetis-dev`, `libhypre-dev`, and `time`, then rebuilt successfully.

## Verification

- Configured and built successfully:
  - `CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release`
  - `cmake --build build -j`
- Smoke-tested the original unit brick path:
  - `mpirun --allow-run-as-root -np 2 ./build/heat_conduction_mpi examples/unit_brick.plc i 6 6 3`
  - Mesh: 196 vertices, 216 wedge cells; HYPRE PCG converged in 6 iterations.
- Validated the new composite geometry/material setup on a reduced mesh:
  - `mpirun --allow-run-as-root -np 2 ./build/heat_conduction_mpi examples/large_composite_heat_multimaterial.plc i 30 25 3`
  - Mesh: 2,880 vertices, 2,850 wedge cells; HYPRE PCG converged in 71 iterations.

## Large benchmark results

Command shape:

```bash
/usr/bin/time -v mpirun --allow-run-as-root -np <ranks> \
  ./build/heat_conduction_mpi examples/large_composite_heat_multimaterial.plc i 110 90 9
```

- Mesh: 9,725 footprint vertices, 11,725 footprint triangles, 97,250 wedge vertices, 105,525 wedge cells.
- Matrix: 942,974 nonzeros.
- HYPRE PCG: 526 iterations.
- Relative L2 error: 0.395812.

| MPI ranks | Internal total (s) | Assembly (s) | Solve (s) | Peak RSS/rank (MB) | External wall time |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 4.871292 | 0.590806 | 3.168944 | 268.699219 | 0:05.23 |
| 2 | 4.725387 | 0.346704 | 3.203277 | 215.355469 | 0:05.07 |
| 4 | 4.570179 | 0.220473 | 3.174343 | 210.203125 | 0:04.94 |

Detailed benchmark logs:

- `logs/benchmarks/large_composite_np1_2026-07-11.log`
- `logs/benchmarks/large_composite_np2_2026-07-11.log`
- `logs/benchmarks/large_composite_np4_2026-07-11.log`
- `logs/benchmarks/large_composite_summary_2026-07-11.md`

## Observations

- Assembly runtime and per-rank memory improve with more MPI ranks after the assembly locality change.
- Overall runtime is still dominated by the current gather-to-rank-0 solve path, so adding ranks primarily improves assembly and memory instead of the HYPRE solve time.

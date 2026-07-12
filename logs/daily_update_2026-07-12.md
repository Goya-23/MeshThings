# Daily update - 2026-07-12

## Large heat-conduction case

- Added `examples/large_composite_heat_multimaterial.plc`, a 20-vertex nonconvex composite footprint extruded to height `2.40`.
- The footprint is represented as one outer boundary assembled from multiple polygonal lobes so it can be processed by the existing PLC reader.
- Added twelve material regions over a default matrix material (`k = 5.0`):
  - copper bar (`k = 385.0`)
  - aluminum fin (`k = 205.0`)
  - ceramic lid (`k = 24.0`)
  - insulation west (`k = 0.08`)
  - graphite sheet (`k = 130.0`)
  - steel mount (`k = 45.0`)
  - silicon tile (`k = 148.0`)
  - polymer gap (`k = 0.22`)
  - brass corner (`k = 109.0`)
  - diamond spreader (`k = 1200.0`)
  - air slot (`k = 0.026`)
  - titanium rib (`k = 22.0`)

## Code quality and performance work

- Limited FEM assembly to wedges touching each rank's owned vertices, avoiding full redundant wedge assembly on every MPI rank.
- Added benchmark instrumentation to `heat_conduction_mpi`:
  - stage timings with max seconds across ranks
  - sparse row and nonzero balance across ranks
  - current RSS and peak RSS high-water marks from `/proc/self/status`
  - a machine-parseable `BENCHMARK` summary line
- Added PLC parser validation for negative hole and material-region counts.

## Verification and benchmark runs

Configured and built successfully:

```bash
CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Smoke-tested the new case:

```bash
mpirun --allow-run-as-root -np 2 ./build/heat_conduction_mpi \
  examples/large_composite_heat_multimaterial.plc i 24 20 3
```

- Mesh: 457 footprint vertices, 557 footprint triangles, 1,828 wedge vertices, 1,671 wedge cells.
- All twelve material regions had nonzero cell counts.
- HYPRE PCG converged in 40 iterations with residual `8.46273e-11`.

Benchmarked the full target mesh:

```bash
mpirun --allow-run-as-root -np <ranks> ./build/heat_conduction_mpi \
  examples/large_composite_heat_multimaterial.plc i 110 90 10
```

- Mesh: 9,721 footprint vertices, 10,931 footprint triangles, 106,931 wedge vertices, 109,310 wedge cells.
- Matrix: 106,931 rows, 1,010,640 nonzeros.
- HYPRE PCG converged in 792 iterations for each benchmark run.
- Relative L2 error vs manufactured solution: `0.368548`.

| MPI ranks | Total runtime (s) | Assembly (s) | Solve (s) | Peak RSS max rank (MiB) | Rows/rank min-max | NNZ/rank min-max |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 6.761749 | 0.545334 | 5.092048 | 283.152344 | 106,931-106,931 | 1,010,640-1,010,640 |
| 2 | 10.533068 | 0.756862 | 7.405666 | 226.437500 | 53,464-53,467 | 498,858-511,782 |
| 4 | 8.163529 | 0.394624 | 5.043006 | 220.855469 | 26,706-26,758 | 230,594-281,423 |

Raw benchmark logs:

- `logs/benchmarks/large_composite_smoke_np2_2026-07-12.log`
- `logs/benchmarks/large_composite_np1_nz10_2026-07-12.log`
- `logs/benchmarks/large_composite_np2_nz10_2026-07-12.log`
- `logs/benchmarks/large_composite_np4_nz10_2026-07-12.log`

## Notes and bottlenecks

- The rank-local assembly improvement reduces redundant assembly work.
- End-to-end runtime does not improve monotonically with process count because mesh generation/OpenVolumeMesh conversion are still replicated on every rank and the HYPRE solve path gathers the sparse system to rank 0.
- Peak RSS per rank decreased from `283.15 MiB` on one rank to `220.86 MiB` on four ranks, but total process memory increases because each rank still keeps a full copy of the mesh.

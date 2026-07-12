# Large composite multi-material heat-conduction benchmark

Case file: `examples/large_composite_heat_multimaterial.plc`

This case is intended to exercise the heat-conduction pipeline on a larger composite solid:

- 20-vertex nonconvex footprint assembled from multiple polygonal lobes and extruded to height `2.40`
- default matrix conductivity `k = 5.0`
- 12 z-aware material regions assigned by wedge centroid
- benchmark target command uses `nx = 110`, `ny = 90`, `nz = 10`, producing `109,310` wedge cells

Material regions are evaluated in file order; later regions override earlier regions for wedges whose centroids lie in overlapping boxes.

## Suggested run

```bash
CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
mpirun --allow-run-as-root -np 4 ./build/heat_conduction_mpi \
  examples/large_composite_heat_multimaterial.plc i 110 90 10
```

The case intentionally omits the command-line conductivity argument so the PLC default conductivity and material regions are used.

## Benchmark record

The benchmark driver prints:

- mesh size and material cell counts
- sparse matrix row and nonzero balance across ranks
- stage timings for parse, triangulation, extrusion, OpenVolumeMesh conversion, partitioning, reordering, assembly, solve, error computation, and total runtime
- maximum current RSS and peak RSS across MPI ranks
- a machine-parseable `BENCHMARK` line

Measured results are recorded in `logs/daily_update_2026-07-12.md` and raw command logs under `logs/benchmarks/`.

## Verification record

Run on 2026-07-12 with command:

```text
mpirun --allow-run-as-root -np <ranks> ./build/heat_conduction_mpi examples/large_composite_heat_multimaterial.plc i 110 90 10
```

Common mesh and solve properties across the benchmark runs:

- Footprint CDT: 9,721 vertices, 10,931 triangles
- Wedge mesh: 106,931 vertices, 109,310 wedge cells
- Matrix: 106,931 rows, 1,010,640 nonzeros
- HYPRE PCG: 792 iterations, final relative residual about `9.1e-11`
- Relative L2 error vs manufactured solution: `0.368548`

| MPI ranks | Total runtime (s) | Assembly (s) | Solve (s) | Peak RSS max rank (MiB) | Rows/rank min-max | NNZ/rank min-max |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 6.761749 | 0.545334 | 5.092048 | 283.152344 | 106,931-106,931 | 1,010,640-1,010,640 |
| 2 | 10.533068 | 0.756862 | 7.405666 | 226.437500 | 53,464-53,467 | 498,858-511,782 |
| 4 | 8.163529 | 0.394624 | 5.043006 | 220.855469 | 26,706-26,758 | 230,594-281,423 |

The benchmark currently remains limited by full mesh construction on each rank and the gather-to-rank-0 HYPRE solve path. The rank-local assembly change reduces redundant assembly work, but end-to-end scaling still depends on replacing the serial solve and replicated preprocessing.

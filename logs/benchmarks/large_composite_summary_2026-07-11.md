# Large composite heat benchmark summary - 2026-07-11

Case: `examples/large_composite_heat_multimaterial.plc`

Command shape:

```bash
/usr/bin/time -v mpirun --allow-run-as-root -np <ranks> \
  ./build/heat_conduction_mpi examples/large_composite_heat_multimaterial.plc i 110 90 9
```

## Mesh and solver

- Footprint CDT: 9,725 vertices, 11,725 triangles.
- Wedge mesh: 97,250 vertices, 105,525 wedge cells.
- OpenVolumeMesh: 105,525 cells, 307,663 faces, 97,250 vertices.
- Sparse matrix: 942,974 nonzeros.
- HYPRE PCG: 526 iterations.
- Relative L2 error vs manufactured solution: 0.395812.

## Process-count comparison

| MPI ranks | Internal total (s) | Assembly (s) | Solve (s) | Peak RSS/rank (MB) | External max RSS (KB) | External wall time |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 4.871292 | 0.590806 | 3.168944 | 268.699219 | 275148 | 0:05.23 |
| 2 | 4.725387 | 0.346704 | 3.203277 | 215.355469 | 220340 | 0:05.07 |
| 4 | 4.570179 | 0.220473 | 3.174343 | 210.203125 | 214820 | 0:04.94 |

## Notes

- The `110 90 9` run satisfies the large-case target with 105,525 wedge elements.
- Assembly time decreased from 0.590806 s at 1 rank to 0.220473 s at 4 ranks.
- Peak per-rank RSS decreased from 268.699219 MB at 1 rank to 210.203125 MB at 4 ranks.
- Solve time stayed near 3.2 s because the current solver implementation gathers the system to rank 0 and solves through the rank-0 HYPRE path.

## Detailed logs

- `logs/benchmarks/large_composite_np1_2026-07-11.log`
- `logs/benchmarks/large_composite_np2_2026-07-11.log`
- `logs/benchmarks/large_composite_np4_2026-07-11.log`

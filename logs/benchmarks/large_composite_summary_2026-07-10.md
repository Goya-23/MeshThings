# Large composite heat benchmark - 2026-07-10

Case: `examples/large_composite_heat_multimaterial.plc`

Command shape:

```bash
/usr/bin/time -v mpirun --allow-run-as-root -np <ranks> \
  ./build/heat_conduction_mpi examples/large_composite_heat_multimaterial.plc i 110 90 9
```

Mesh:

- Footprint: 9,725 vertices, 11,725 triangles
- Wedge mesh: 97,250 vertices, 105,525 wedge cells
- OpenVolumeMesh: 105,525 cells, 307,663 faces, 97,250 vertices
- Sparse matrix nonzeros: 942,974
- HYPRE PCG iterations: 526
- Relative L2 error vs manufactured solution: 0.395812

| MPI ranks | Internal total (s) | Assembly (s) | Solve (s) | Internal peak RSS/rank (MB) | `/usr/bin/time` wall | `/usr/bin/time` max RSS (KB) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 4.757378 | 0.519714 | 3.135757 | 268.804688 | 0:05.13 | 275256 |
| 2 | 4.570183 | 0.285066 | 3.122220 | 215.320312 | 0:04.93 | 220188 |
| 4 | 4.568784 | 0.183612 | 3.158071 | 210.800781 | 0:04.95 | 215396 |

Notes:

- Assembly time decreases from 0.519714 s on 1 rank to 0.183612 s on 4 ranks after assembling only wedges incident to each rank's owned vertices.
- Peak per-rank RSS drops from 268.804688 MB on 1 rank to 210.800781 MB on 4 ranks.
- Total runtime is nearly flat at 2-4 ranks because `HypreSolver::solve` still gathers all rows to rank 0 and solves with `MPI_COMM_SELF`.

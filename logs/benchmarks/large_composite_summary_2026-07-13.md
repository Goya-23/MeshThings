# Large composite heat benchmark - 2026-07-13

Command template:

```bash
mpirun --allow-run-as-root -np <N> ./build/heat_conduction_mpi \
  examples/large_composite_heat_multimaterial.plc i 110 90 10
```

Mesh:

- Footprint CDT: 9,721 vertices, 10,623 triangles
- Wedge mesh: 106,931 vertices, 106,230 wedge cells
- Matrix: 106,931 rows, 992,274 nonzeros
- Relative L2 error: 0.355407

| MPI processes | Total seconds | Peak RSS max (MB) | `/usr/bin/time` elapsed | `/usr/bin/time` max RSS (KB) | Raw log |
| ---: | ---: | ---: | ---: | ---: | --- |
| 1 | 4.95674 | 277.492 | 0:05.30 | 284048 | `large_composite_np1_2026-07-13.log` |
| 2 | 4.80963 | 222.293 | 0:05.15 | 227560 | `large_composite_np2_2026-07-13.log` |
| 4 | 4.74083 | 218.25 | 0:05.10 | 223152 | `large_composite_np4_2026-07-13.log` |

The assembly stage scales down from 0.511282 seconds at 1 rank to 0.183385 seconds at 4 ranks. End-to-end runtime changes less because the current solve path gathers the matrix to rank 0 before invoking HYPRE with `MPI_COMM_SELF`.

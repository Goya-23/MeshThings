# Large composite heat benchmark - 2026-07-15

Case: `examples/large_composite_heat_multimaterial.plc`

Mesh arguments: `i 110 90 10`

Raw output files:

- `logs/benchmarks/large_composite_heat_multimaterial_2026-07-15_np1.txt`
- `logs/benchmarks/large_composite_heat_multimaterial_2026-07-15_np2.txt`
- `logs/benchmarks/large_composite_heat_multimaterial_2026-07-15_np4.txt`

## Mesh and material summary

All runs used:

```text
mpirun --allow-run-as-root -np <N> ./build/heat_conduction_mpi examples/large_composite_heat_multimaterial.plc i 110 90 10
```

Shared mesh statistics:

- Footprint CDT: 9,723 vertices, 11,065 triangles.
- Wedge mesh: 106,953 vertices, 110,650 wedge cells.
- OpenVolumeMesh: 110,650 cells, 322,065 faces, 106,953 vertices.
- Relative L2 error: `0.385714`.

Material cell counts:

| Material | Wedge cells |
| --- | ---: |
| default | 21597 |
| copper_bus_left | 6922 |
| copper_bus_right | 3508 |
| aluminum_frame | 11660 |
| ceramic_cap_upper | 16374 |
| graphite_midplane | 4568 |
| insulation_notch | 11135 |
| insulation_bridge | 15758 |
| steel_fastener_a | 2450 |
| steel_fastener_b | 690 |
| diamond_spreader | 6476 |
| polymer_skin_top | 7104 |
| coolant_channel | 2408 |

## Runtime and memory comparison

| MPI processes | Total runtime (s) | CDT (s) | Extrude/OVM (s) | Partition (s) | Reorder (s) | Assembly (s) | Solve (s) | Gather/error (s) | Max peak RSS (KB) | `/usr/bin/time` max RSS (KB) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 4.46046 | 0.51864 | 0.583542 | 0.000835059 | 0.00517868 | 0.524808 | 2.81589 | 0.0107197 | 291704 | 291704 |
| 2 | 4.29625 | 0.524687 | 0.601014 | 0.0572942 | 0.00585242 | 0.276725 | 2.81518 | 0.0173047 | 232700 | 232700 |
| 4 | 4.44741 | 0.530986 | 0.667439 | 0.113618 | 0.0097825 | 0.213902 | 2.93304 | 0.058876 | 229748 | 229748 |

| MPI processes | Rows min/max | Sparse nnz min/max | Incident wedges assembled min/max | HYPRE PCG iterations |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 106953 / 106953 | 1020116 / 1020116 | 110650 / 110650 | 427 |
| 2 | 53474 / 53479 | 507870 / 512246 | 54290 / 56380 | 427 |
| 4 | 24833 / 28614 | 186765 / 312443 | 19130 / 35051 | 428 |

## Notes

- The FEM assembly improvement reduces per-rank incident wedge assembly from all 110,650 cells on each rank to the partition-local incident subset.
- Peak RSS drops from 291,704 KB with one process to 229,748 KB with four processes.
- End-to-end runtime remains limited by replicated input/CDT/extrusion work and the current gather-to-rank-0 HYPRE solve path.

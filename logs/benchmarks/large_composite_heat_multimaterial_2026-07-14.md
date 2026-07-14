# Large composite heat-conduction benchmark - 2026-07-14

Case: `examples/large_composite_heat_multimaterial.plc`

Resolution: `nx = 110`, `ny = 90`, `nz = 10`

## Commands

```bash
CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

for np in 1 2 4; do
  /usr/bin/time -v mpirun --allow-run-as-root -np "${np}" \
    ./build/heat_conduction_mpi examples/large_composite_heat_multimaterial.plc i 110 90 10
done
```

## Results

Build:

```text
CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Small regression:

```text
mpirun --allow-run-as-root -np 2 ./build/heat_conduction_mpi examples/unit_brick.plc i 6 6 3
```

The regression completed successfully with HYPRE PCG converging in 6 iterations
and the new `BENCHMARK` output emitted.

Large-case mesh:

- footprint CDT: 9,723 vertices, 11,065 triangles
- wedge mesh: 106,953 vertices, 110,650 wedge cells
- OpenVolumeMesh: 110,650 cells, 322,065 faces, 106,953 vertices

Material cell counts:

| Material | Conductivity | Wedge cells |
| --- | ---: | ---: |
| default | 7.5 | 21,597 |
| copper_bus_left | 385 | 6,922 |
| copper_bus_right | 385 | 3,508 |
| aluminum_frame | 205 | 11,660 |
| ceramic_cap_upper | 24 | 16,374 |
| graphite_midplane | 130 | 4,568 |
| insulation_notch | 0.06 | 11,135 |
| insulation_bridge | 0.08 | 15,758 |
| steel_fastener_a | 16 | 2,450 |
| steel_fastener_b | 16 | 690 |
| diamond_spreader | 900 | 6,476 |
| polymer_skin_top | 0.22 | 7,104 |
| coolant_channel | 0.60 | 2,408 |

Process-count comparison:

| MPI processes | Driver total (s) | `/usr/bin/time` wall (s) | Assembly (s) | Solve (s) | Peak RSS max (KB) | Rows min/max | NNZ min/max | Incident wedges min/max |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 4.7185 | 5.12 | 0.616573 | 2.8574 | 291728 | 106953 / 106953 | 1020116 / 1020116 | 110650 / 110650 |
| 2 | 4.5285 | 4.92 | 0.369824 | 2.83551 | 232888 | 53474 / 53479 | 507870 / 512246 | 54290 / 56380 |
| 4 | 4.92442 | 5.35 | 0.294131 | 3.25617 | 229872 | 24833 / 28614 | 186765 / 312443 | 19130 / 35051 |

Raw `BENCHMARK` lines:

```text
BENCHMARK processes=1 nx=110 ny=90 nz=10 vertices=106953 wedges=110650 input_cdt_s=0.53939 extrude_ovm_s=0.682804 partition_s=0.00125271 reorder_s=0.00903991 assembly_s=0.616573 solve_s=2.8574 gather_l2_error_s=0.0110146 total_s=4.7185 rows_min=106953 rows_max=106953 nnz_min=1020116 nnz_max=1020116 incident_wedges_min=110650 incident_wedges_max=110650 rss_peak_kb_max=291728 l2_error=0.385714
BENCHMARK processes=2 nx=110 ny=90 nz=10 vertices=106953 wedges=110650 input_cdt_s=0.535215 extrude_ovm_s=0.693866 partition_s=0.0621603 reorder_s=0.00779359 assembly_s=0.369824 solve_s=2.83551 gather_l2_error_s=0.0199546 total_s=4.5285 rows_min=53474 rows_max=53479 nnz_min=507870 nnz_max=512246 incident_wedges_min=54290 incident_wedges_max=56380 rss_peak_kb_max=232888 l2_error=0.385714
BENCHMARK processes=4 nx=110 ny=90 nz=10 vertices=106953 wedges=110650 input_cdt_s=0.566923 extrude_ovm_s=0.800674 partition_s=0.173994 reorder_s=0.0157295 assembly_s=0.294131 solve_s=3.25617 gather_l2_error_s=0.0284993 total_s=4.92442 rows_min=24833 rows_max=28614 nnz_min=186765 nnz_max=312443 incident_wedges_min=19130 incident_wedges_max=35051 rss_peak_kb_max=229872 l2_error=0.385714
```

Raw command logs:

- `logs/benchmarks/large_composite_heat_multimaterial_2026-07-14_np1.txt`
- `logs/benchmarks/large_composite_heat_multimaterial_2026-07-14_np2.txt`
- `logs/benchmarks/large_composite_heat_multimaterial_2026-07-14_np4.txt`

## Notes

- The case satisfies the requested size target with 110,650 wedge cells.
- Assembly time decreases as the process count increases because each rank now
  assembles only wedges incident to owned vertices.
- Total runtime remains dominated by duplicated mesh construction on every rank
  and by the current serial HYPRE solve after gathering the matrix to rank 0.
- The L2 error is reported for continuity with the existing manufactured
  solution check; because the case uses strong material contrasts, it is a
  smoke-test metric rather than a strict multi-material accuracy proof.

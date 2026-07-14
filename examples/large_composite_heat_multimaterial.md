# Large composite multi-material heat-conduction case

Case file: `examples/large_composite_heat_multimaterial.plc`

This case is intended for stress-testing the heat-conduction pipeline above
100,000 wedge elements. The 2D footprint is a stitched, non-rectangular
polygonal silhouette that represents a composite made from several arbitrary
polygonal lobes. It is extruded through a thickness of `2.40` in `z`.

## Geometry and materials

- 22-vertex concave outer footprint
- extrusion height `H = 2.40`
- default matrix conductivity `k = 7.5`
- 12 centroid-assigned material regions:
  - `copper_bus_left`, `k = 385.0`
  - `copper_bus_right`, `k = 385.0`
  - `aluminum_frame`, `k = 205.0`
  - `ceramic_cap_upper`, `k = 24.0`
  - `graphite_midplane`, `k = 130.0`
  - `insulation_notch`, `k = 0.06`
  - `insulation_bridge`, `k = 0.08`
  - `steel_fastener_a`, `k = 16.0`
  - `steel_fastener_b`, `k = 16.0`
  - `diamond_spreader`, `k = 900.0`
  - `polymer_skin_top`, `k = 0.22`
  - `coolant_channel`, `k = 0.60`

Later material regions in the PLC override earlier ones for wedges whose
centroids lie in more than one material box.

## Benchmark command

```bash
CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

for np in 1 2 4; do
  /usr/bin/time -v mpirun --allow-run-as-root -np "${np}" \
    ./build/heat_conduction_mpi examples/large_composite_heat_multimaterial.plc i 110 90 10
done
```

The command intentionally omits the final conductivity override argument so the
PLC default and material-region conductivities are used.

## Verification and benchmark record

Run on 2026-07-14 with `nx = 110`, `ny = 90`, `nz = 10`:

- footprint CDT: 9,723 vertices, 11,065 triangles
- wedge mesh: 106,953 vertices, 110,650 wedge cells
- HYPRE PCG: 427-428 iterations
- relative L2 error vs the manufactured solution: `0.385714`

| MPI processes | Driver total (s) | `/usr/bin/time` wall (s) | Assembly (s) | Solve (s) | Peak RSS max (KB) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 4.7185 | 5.12 | 0.616573 | 2.8574 | 291728 |
| 2 | 4.5285 | 4.92 | 0.369824 | 2.83551 | 232888 |
| 4 | 4.92442 | 5.35 | 0.294131 | 3.25617 | 229872 |

The per-rank assembly time and peak resident set size improve with additional
processes. End-to-end runtime does not scale monotonically because each rank
still builds the full mesh and the current HYPRE path gathers the matrix for a
serial solve on rank 0.

Detailed logs:

- `logs/benchmarks/large_composite_heat_multimaterial_2026-07-14.md`
- `logs/benchmarks/large_composite_heat_multimaterial_2026-07-14_np1.txt`
- `logs/benchmarks/large_composite_heat_multimaterial_2026-07-14_np2.txt`
- `logs/benchmarks/large_composite_heat_multimaterial_2026-07-14_np4.txt`

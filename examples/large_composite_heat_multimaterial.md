# Large composite multi-material heat-conduction case

Case file: `examples/large_composite_heat_multimaterial.plc`

This case is intended for large heat-conduction performance runs:

- 24-vertex concave footprint assembled from several polygonal protrusions and notches.
- Extrusion height `2.40`, using PLC-provided z thickness.
- Default composite matrix conductivity `k = 2.5`.
- Twelve overlapping material boxes assigned by wedge centroid:
  - `copper_bus`, `k = 385.0`
  - `aluminum_frame`, `k = 205.0`
  - `ceramic_plate`, `k = 24.0`
  - `graphite_spreader`, `k = 130.0`
  - `insulation_slot`, `k = 0.06`
  - `steel_anchor`, `k = 45.0`
  - `silicon_die`, `k = 148.0`
  - `polymer_lid`, `k = 0.22`
  - `coolant_channel`, `k = 0.60`
  - `thermal_interface`, `k = 8.0`
  - `tungsten_insert`, `k = 170.0`
  - `air_gap`, `k = 0.026`

Later material regions in the PLC override earlier regions when a wedge centroid lies in more than one material box.

## Suggested large run

```bash
CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
mpirun --allow-run-as-root -np 4 ./build/heat_conduction_mpi \
  examples/large_composite_heat_multimaterial.plc i 110 90 9
```

The case intentionally omits the command-line conductivity argument so the PLC default conductivity and material regions are used.

## Benchmark record

Run on 2026-07-11 with `nx=110`, `ny=90`, `nz=9`:

```text
Footprint CDT: 9725 vertices, 11725 triangles
Wedge mesh: 97250 vertices, 105525 wedge cells
OpenVolumeMesh cells: 105525, faces: 307663, vertices: 97250
Sparse matrix nnz: 942974
HYPRE PCG iterations: 526
Relative L2 error vs manufactured solution: 0.395812
```

Process-count comparison:

| MPI ranks | Internal total (s) | Assembly (s) | Solve (s) | Peak RSS/rank (MB) | External wall time |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 4.871292 | 0.590806 | 3.168944 | 268.699219 | 0:05.23 |
| 2 | 4.725387 | 0.346704 | 3.203277 | 215.355469 | 0:05.07 |
| 4 | 4.570179 | 0.220473 | 3.174343 | 210.203125 | 0:04.94 |

Detailed logs:

- `logs/benchmarks/large_composite_np1_2026-07-11.log`
- `logs/benchmarks/large_composite_np2_2026-07-11.log`
- `logs/benchmarks/large_composite_np4_2026-07-11.log`

Assembly runtime and peak per-rank memory improve as ranks increase. End-to-end runtime changes modestly because the current HYPRE path still gathers the assembled system to rank 0 before solving.

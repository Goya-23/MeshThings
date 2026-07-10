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

Run on 2026-07-10 with `nx=110`, `ny=90`, `nz=9`:

```text
Footprint CDT: 9725 vertices, 11725 triangles
Wedge mesh: 97250 vertices, 105525 wedge cells
OpenVolumeMesh cells: 105525, faces: 307663, vertices: 97250
Sparse matrix nnz: 942974
HYPRE PCG iterations: 526
Relative L2 error vs manufactured solution: 0.395812
```

Process-count comparison:

| MPI ranks | Total solver-reported time (s) | Assembly time (s) | Solve time (s) | Peak RSS/rank (MB) | External wall time |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 4.757378 | 0.519714 | 3.135757 | 268.804688 | 0:05.13 |
| 2 | 4.570183 | 0.285066 | 3.122220 | 215.320312 | 0:04.93 |
| 4 | 4.568784 | 0.183612 | 3.158071 | 210.800781 | 0:04.95 |

Detailed logs:

- `logs/benchmarks/large_composite_np1_2026-07-10.log`
- `logs/benchmarks/large_composite_np2_2026-07-10.log`
- `logs/benchmarks/large_composite_np4_2026-07-10.log`

Assembly scales down with additional ranks after localizing element contributions by owned rows. End-to-end runtime remains mostly flat because the current solver path gathers the matrix to rank 0 and runs HYPRE on `MPI_COMM_SELF`.

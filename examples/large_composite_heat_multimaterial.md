# Large composite multi-material heat-conduction benchmark

Case file: `examples/large_composite_heat_multimaterial.plc`

This benchmark is intended to stress the full 3D heat-conduction pipeline with a non-rectangular 2D footprint, z extrusion, material overrides, METIS partitioning, FEM assembly, and the HYPRE solve path.

## Geometry and materials

- 20-vertex concave composite footprint, representing the union outline of several offset polygonal plates.
- Extrusion height: `2.40`.
- Default matrix conductivity: `k = 6.5`.
- Material regions are axis-aligned 3D boxes assigned by wedge centroid:
  - `copper_busbar`, `k = 385.0`
  - `graphite_spreader`, `k = 130.0`
  - `ceramic_bridge`, `k = 24.0`
  - `polymer_slot`, `k = 0.18`
  - `steel_frame`, `k = 45.0`
  - `aluminum_fin`, `k = 205.0`
  - `titanium_mount`, `k = 22.0`
  - `foam_relief`, `k = 0.045`
  - `silicon_die`, `k = 149.0`
  - `air_gap`, `k = 0.026`
  - `silver_via_field`, `k = 429.0`
  - `thermal_grease_layer`, `k = 8.5`

Later material regions in the PLC override earlier regions when centroid boxes overlap.

## Large benchmark command

The verification command intentionally omits the optional uniform conductivity argument so the PLC material map is used:

```bash
mpirun --allow-run-as-root -np 1 ./build/heat_conduction_mpi \
  examples/large_composite_heat_multimaterial.plc i 110 90 10
```

The same command should be run with `-np 2` and `-np 4` for process-count comparisons. The `110 90 10` grid/layer setting produces 106,230 wedge cells on this irregular footprint.

## Benchmark output

`heat_conduction_mpi` prints:

- stage timings as `TIMING <stage> seconds min/avg/max = ...`
- local sparse-matrix row/nnz balance
- process RSS and peak RSS summaries
- a machine-parseable `BENCHMARK` line with process count, mesh size, total wall time, peak RSS, and L2 error

Detailed run records are stored under `logs/benchmarks/`, with a daily summary in `logs/daily_update_2026-07-13.md`.

## Verification record

Run on 2026-07-13 with:

```bash
CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
mpirun --allow-run-as-root -np <N> ./build/heat_conduction_mpi \
  examples/large_composite_heat_multimaterial.plc i 110 90 10
```

Mesh and solve invariants:

- Footprint CDT: 9,721 vertices, 10,623 triangles
- Wedge mesh: 106,931 vertices, 106,230 wedge cells
- Matrix: 106,931 rows, 992,274 nonzeros
- HYPRE PCG: 539 iterations
- Relative L2 error: 0.355407

| MPI processes | Total seconds | Peak RSS max (MB) | Assembly seconds max | HYPRE solve seconds max | Raw log |
| ---: | ---: | ---: | ---: | ---: | --- |
| 1 | 4.960215 | 277.566406 | 0.506088 | 3.380463 | `logs/benchmarks/large_composite_np1_2026-07-13.log` |
| 2 | 4.803107 | 221.656250 | 0.276566 | 3.382213 | `logs/benchmarks/large_composite_np2_2026-07-13.log` |
| 4 | 4.740419 | 218.136719 | 0.179491 | 3.411262 | `logs/benchmarks/large_composite_np4_2026-07-13.log` |

Assembly time decreases as the incident-wedge assembly work is split across ranks. End-to-end runtime is only modestly lower because the current HYPRE path gathers the matrix to rank 0 and solves with `MPI_COMM_SELF`.

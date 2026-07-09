# Daily update - 2026-07-09

## Large composite heat-conduction case

- Added `examples/large_composite_heat_multimaterial.plc`, a 25-vertex nonconvex composite footprint with extrusion height `2.40`.
- The footprint represents joined heat-spreader, busbar, mounting-tab, and cooling-lip polygonal subdomains.
- Added 12 layered or overlapping material regions over a default matrix material:
  - aluminum frame (`k = 205.0`)
  - copper busbars (`k = 385.0`)
  - ceramic power stage (`k = 24.0`)
  - graphite spreader (`k = 130.0`)
  - aerogel slot (`k = 0.035`)
  - polymer gap (`k = 0.18`)
  - steel mount (`k = 45.0`)
  - silicon die (`k = 149.0`)
  - coolant channel (`k = 0.60`)
  - diamond insert (`k = 1000.0`)
  - low-k lid (`k = 0.12`)

## Code improvements

- Added a dense clipped-structured triangulation path for large arbitrary footprints so non-rectangular cases can be generated quickly at high resolution.
- Added topological boundary-edge detection before wedge extrusion so clipped arbitrary footprints get Dirichlet boundary metadata on the generated mesh boundary.
- Reduced FEM assembly work per MPI rank by assembling only wedges that touch vertices owned by that rank.
- Added `METRIC` output lines to the MPI demo for phase runtime and max/sum resident set size across ranks.
- Added `scripts/benchmark_heat_case.sh` for repeatable process-count comparisons.

## Verification

- Installed documented system dependencies in this environment:
  - `openmpi-bin libopenmpi-dev libmetis-dev libhypre-dev`
- Configured and built successfully:
  - `CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release`
  - `cmake --build build -j`
- Re-ran compatibility cases:
  - `mpirun --allow-run-as-root -np 2 ./build/heat_conduction_mpi examples/unit_brick.plc i 6 6 3`
  - `mpirun --allow-run-as-root -np 2 ./build/heat_conduction_mpi examples/complex_heat_multimaterial.plc i 18 18 6`

## Large-case process-count comparison

Command:

```bash
scripts/benchmark_heat_case.sh "1 2 4" examples/large_composite_heat_multimaterial.plc 110 90 8
```

Mesh:

- Footprint: 7069 vertices, 13688 triangles.
- Wedge mesh: 63621 vertices, 109504 wedge cells.
- OpenVolumeMesh: 109504 cells, 289240 faces, 63621 vertices.
- HYPRE PCG: 10 iterations.
- Relative L2 error: `0.310276`.

Summary:

| MPI processes | total seconds max | max RSS KiB | summed RSS KiB |
| ---: | ---: | ---: | ---: |
| 1 | 1.424189 | 234160 | 234160 |
| 2 | 1.294379 | 187488 | 352356 |
| 4 | 1.129636 | 182988 | 603796 |

Full benchmark logs:

- `logs/benchmarks/heat_np1_nx110_ny90_nz8.log`
- `logs/benchmarks/heat_np2_nx110_ny90_nz8.log`
- `logs/benchmarks/heat_np4_nx110_ny90_nz8.log`
- `logs/benchmarks/large_composite_summary_2026-07-09.csv`

# Daily update - 2026-07-15

## Heat-conduction case

- Added `examples/large_composite_heat_multimaterial.plc`, a large 22-vertex stitched polygonal footprint with extrusion height `2.40`.
- The case combines arbitrary 2D polygonal lobes into one extruded footprint and assigns 12 material regions over a default matrix material.
- Benchmark mesh target: `nx = 110`, `ny = 90`, `nz = 10`, expected to produce more than 100,000 wedge elements.

## Code improvements

- Added per-stage benchmark timing to `heat_conduction_mpi`.
- Added per-rank owned-row, sparse-nnz, incident-wedge, and peak-RSS balance reporting.
- Improved FEM assembly so each rank assembles only wedges incident to its owned vertices instead of redundantly assembling the full mesh on every rank.

## Verification

- Installed missing runtime/build packages in the cloud image: `openmpi-bin`, `libopenmpi-dev`, `libmetis-dev`, `libhypre-dev`, and `time`.
- Configured and built successfully:
  - `CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release`
  - `cmake --build build -j`
- Re-ran the original unit brick path for backward compatibility:
  - `mpirun --allow-run-as-root -np 2 ./build/heat_conduction_mpi examples/unit_brick.plc i 6 6 3`
  - Mesh: 196 wedge vertices, 216 wedge cells.
  - HYPRE PCG converged in 6 iterations with residual `1.72381e-12`.
- Ran the large composite case with `nx = 110`, `ny = 90`, `nz = 10`:
  - Mesh: 9,723 footprint vertices, 11,065 footprint triangles, 106,953 wedge vertices, 110,650 wedge cells.
  - Relative L2 error was `0.385714` for 1, 2, and 4 MPI processes.
  - Runtime / max peak RSS comparison:
    - 1 process: `4.46046 s`, `291704 KB`.
    - 2 processes: `4.29625 s`, `232700 KB`.
    - 4 processes: `4.44741 s`, `229748 KB`.
  - Raw outputs and detailed timing breakdowns are recorded in `logs/benchmarks/large_composite_heat_multimaterial_2026-07-15*.txt` and `logs/benchmarks/large_composite_heat_multimaterial_2026-07-15.md`.

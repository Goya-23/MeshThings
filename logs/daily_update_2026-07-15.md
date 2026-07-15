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

- Benchmark results will be recorded after build and MPI runs complete.

# Daily update - 2026-07-13

## Large heat-conduction benchmark case

- Added `examples/large_composite_heat_multimaterial.plc`, a 20-vertex concave composite footprint extruded to height `2.40`.
- Added twelve material regions over a default matrix material, including copper, graphite, ceramic, polymer, steel, aluminum, titanium, foam, silicon, air, silver, and grease regions.
- Target verification command:
  - `mpirun --allow-run-as-root -np <N> ./build/heat_conduction_mpi examples/large_composite_heat_multimaterial.plc i 110 90 8`
- The chosen mesh density is intended to exceed 100,000 wedge cells.

## Code quality and benchmark instrumentation

- Reduced FEM assembly work by assembling only wedges incident to vertices owned by the current MPI rank.
- Added benchmark telemetry to `examples/heat_conduction_mpi.cpp`:
  - per-stage min/average/max timings across ranks
  - local sparse matrix row and nonzero balance
  - RSS and peak RSS memory summaries
  - machine-parseable `BENCHMARK` output

## Verification and process-count comparison

Benchmark runs will be recorded after the pre-test implementation commit.

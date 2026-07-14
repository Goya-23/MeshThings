# Daily update - 2026-07-14

## Heat-conduction case

- Added `examples/large_composite_heat_multimaterial.plc`, a 22-vertex
  non-rectangular composite footprint extruded to height `2.40`.
- Added 12 material regions over a default matrix material, covering high
  conductivity copper/diamond spreaders, moderate metal/ceramic/graphite
  regions, and low conductivity insulation/polymer/coolant regions.
- Target benchmark resolution is `nx = 110`, `ny = 90`, `nz = 10`, which is
  intended to exercise more than 100,000 wedge cells.

## Code quality and instrumentation

- Updated MPI heat-driver output with stage timings, row/nnz/incident-wedge
  balance, peak RSS memory summaries, and a machine-parseable `BENCHMARK` line.
- Scoped FEM assembly to wedges incident to each rank's owned vertices instead
  of assembling the full global sparse map on every rank.

## Verification

- Build, run, and process-count benchmark results will be recorded in
  `logs/benchmarks/large_composite_heat_multimaterial_2026-07-14.md`.

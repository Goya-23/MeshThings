# Daily update - 2026-07-12

## Large heat-conduction case

- Added `examples/large_composite_heat_multimaterial.plc`, a 20-vertex nonconvex composite footprint extruded to height `2.40`.
- The footprint is represented as one outer boundary assembled from multiple polygonal lobes so it can be processed by the existing PLC reader.
- Added twelve material regions over a default matrix material (`k = 5.0`):
  - copper bar (`k = 385.0`)
  - aluminum fin (`k = 205.0`)
  - ceramic lid (`k = 24.0`)
  - insulation west (`k = 0.08`)
  - graphite sheet (`k = 130.0`)
  - steel mount (`k = 45.0`)
  - silicon tile (`k = 148.0`)
  - polymer gap (`k = 0.22`)
  - brass corner (`k = 109.0`)
  - diamond spreader (`k = 1200.0`)
  - air slot (`k = 0.026`)
  - titanium rib (`k = 22.0`)

## Code quality and performance work

- Limited FEM assembly to wedges touching each rank's owned vertices, avoiding full redundant wedge assembly on every MPI rank.
- Added benchmark instrumentation to `heat_conduction_mpi`:
  - stage timings with max seconds across ranks
  - sparse row and nonzero balance across ranks
  - current RSS and peak RSS high-water marks from `/proc/self/status`
  - a machine-parseable `BENCHMARK` summary line
- Added PLC parser validation for negative hole and material-region counts.

## Verification and benchmark runs

Pending first benchmark run after the implementation commit.

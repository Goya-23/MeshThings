# Large composite multi-material heat-conduction case

Case file: `examples/large_composite_heat_multimaterial.plc`

This case exercises a denser heat-conduction workload than the earlier complex example:

- 25-vertex nonconvex footprint representing joined arbitrary 2D polygonal subdomains.
- z extrusion height `2.40`.
- default matrix conductivity `k = 2.5`.
- 12 layered or overlapping material boxes assigned by wedge centroid:
  - aluminum frame, copper busbars, ceramic power stage, graphite spreader, aerogel slot,
    polymer gap, steel mount, silicon die, coolant channel, diamond insert, and low-k lid.
- Later material boxes in the PLC override earlier boxes in overlap regions.

## Suggested large run

```bash
CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
mpirun --allow-run-as-root -np 4 ./build/heat_conduction_mpi \
  examples/large_composite_heat_multimaterial.plc i 110 90 8
```

The `110 x 90 x 8` settings intentionally produce more than 100,000 wedge cells for this
composite footprint. The executable reports `METRIC` lines with per-phase max/average
runtime and max/sum resident set size across MPI ranks.

## Benchmark helper

Use `scripts/benchmark_heat_case.sh` to run the same case with multiple process counts:

```bash
scripts/benchmark_heat_case.sh "1 2 4" examples/large_composite_heat_multimaterial.plc 110 90 8
```

The helper leaves the complete command output in `logs/benchmarks/` and prints a compact
summary containing process count, wedge cells, total runtime, and memory.

## Verification and process-count comparison

Run on 2026-07-09:

```bash
scripts/benchmark_heat_case.sh "1 2 4" examples/large_composite_heat_multimaterial.plc 110 90 8
```

Common mesh and solve results:

```text
Footprint CDT: 7069 vertices, 13688 triangles
Wedge mesh: 63621 vertices, 109504 wedge cells
OpenVolumeMesh cells: 109504, faces: 289240, vertices: 63621
HYPRE PCG iterations: 10
Relative L2 error vs manufactured solution: 0.310276
```

Material cell counts:

| material | conductivity | wedge cells |
| --- | ---: | ---: |
| default | 2.5 | 18196 |
| aluminum_frame | 205.0 | 32289 |
| copper_bus_left | 385.0 | 7496 |
| copper_bus_right | 385.0 | 10408 |
| ceramic_power_stage | 24.0 | 11848 |
| graphite_spreader | 130.0 | 3375 |
| aerogel_slot | 0.035 | 5884 |
| polymer_gap | 0.18 | 4056 |
| steel_mount | 45.0 | 6607 |
| silicon_die | 149.0 | 1611 |
| coolant_channel | 0.60 | 3390 |
| diamond_insert | 1000.0 | 1520 |
| low_k_lid | 0.12 | 2824 |

Process-count comparison from the executable's `METRIC phase=total` line:

| MPI processes | wedge cells | total seconds max | max RSS KiB | summed RSS KiB | log |
| ---: | ---: | ---: | ---: | ---: | --- |
| 1 | 109504 | 1.424189 | 234160 | 234160 | `logs/benchmarks/heat_np1_nx110_ny90_nz8.log` |
| 2 | 109504 | 1.294379 | 187488 | 352356 | `logs/benchmarks/heat_np2_nx110_ny90_nz8.log` |
| 4 | 109504 | 1.129636 | 182988 | 603796 | `logs/benchmarks/heat_np4_nx110_ny90_nz8.log` |

The assembly phase benefits most from more ranks after the per-rank assembly change
(`0.584159s` at 1 rank, `0.276112s` at 2 ranks, `0.140082s` at 4 ranks). The solve still
gathers to rank 0 before using HYPRE, so total runtime does not scale linearly and summed
RSS increases with process count.

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

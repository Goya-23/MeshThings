# Large composite multi-material heat-conduction benchmark

Case file: `examples/large_composite_heat_multimaterial.plc`

This case is intended to exercise the heat-conduction pipeline on a larger composite solid:

- 20-vertex nonconvex footprint assembled from multiple polygonal lobes and extruded to height `2.40`
- default matrix conductivity `k = 5.0`
- 12 z-aware material regions assigned by wedge centroid
- benchmark target command uses `nx = 110`, `ny = 90`, `nz = 9`, producing more than 100,000 wedge cells

Material regions are evaluated in file order; later regions override earlier regions for wedges whose centroids lie in overlapping boxes.

## Suggested run

```bash
CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
mpirun --allow-run-as-root -np 4 ./build/heat_conduction_mpi \
  examples/large_composite_heat_multimaterial.plc i 110 90 9
```

The case intentionally omits the command-line conductivity argument so the PLC default conductivity and material regions are used.

## Benchmark record

The benchmark driver prints:

- mesh size and material cell counts
- sparse matrix row and nonzero balance across ranks
- stage timings for parse, triangulation, extrusion, OpenVolumeMesh conversion, partitioning, reordering, assembly, solve, error computation, and total runtime
- maximum current RSS and peak RSS across MPI ranks
- a machine-parseable `BENCHMARK` line

Measured results are recorded in `logs/daily_update_2026-07-12.md` and raw command logs under `logs/benchmarks/`.

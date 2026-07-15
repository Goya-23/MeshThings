# Large composite multi-material heat-conduction case

Case file: `examples/large_composite_heat_multimaterial.plc`

This case is intended as a larger performance and memory benchmark for the wedge heat-conduction pipeline:

- 22-vertex stitched polygonal footprint that represents a union-like combination of arbitrary 2D lobes.
- Extrusion height `2.40` with z-dependent material boxes.
- Default matrix conductivity `k = 7.5`.
- 12 material regions with overlapping boxes; later records override earlier records by wedge centroid.
- Benchmark resolution: `nx = 110`, `ny = 90`, `nz = 10`, which is expected to produce more than 100,000 wedge elements.

## Material regions

| Material | Conductivity |
| --- | ---: |
| `copper_bus_left` | 385.0 |
| `copper_bus_right` | 385.0 |
| `aluminum_frame` | 205.0 |
| `ceramic_cap_upper` | 24.0 |
| `graphite_midplane` | 130.0 |
| `insulation_notch` | 0.06 |
| `insulation_bridge` | 0.08 |
| `steel_fastener_a` | 16.0 |
| `steel_fastener_b` | 16.0 |
| `diamond_spreader` | 900.0 |
| `polymer_skin_top` | 0.22 |
| `coolant_channel` | 0.60 |

## Benchmark command

```bash
CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
mpirun --allow-run-as-root -np 4 ./build/heat_conduction_mpi \
  examples/large_composite_heat_multimaterial.plc i 110 90 10
```

Use 1, 2, and 4 MPI processes with the same mesh arguments to compare runtime and memory usage. Raw benchmark output for the 2026-07-15 run is recorded under `logs/benchmarks/`.

## Verification record

Run on 2026-07-15:

```text
Results will be populated after the benchmark suite completes.
```

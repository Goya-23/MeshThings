# Large composite multi-material heat-conduction case

Case file: `examples/large_composite_heat_multimaterial.plc`

This case is intended for large heat-conduction performance runs:

- 24-vertex concave footprint assembled from several polygonal protrusions and notches.
- Extrusion height `2.40`, using PLC-provided z thickness.
- Default composite matrix conductivity `k = 2.5`.
- Twelve overlapping material boxes assigned by wedge centroid:
  - `copper_bus`, `k = 385.0`
  - `aluminum_frame`, `k = 205.0`
  - `ceramic_plate`, `k = 24.0`
  - `graphite_spreader`, `k = 130.0`
  - `insulation_slot`, `k = 0.06`
  - `steel_anchor`, `k = 45.0`
  - `silicon_die`, `k = 148.0`
  - `polymer_lid`, `k = 0.22`
  - `coolant_channel`, `k = 0.60`
  - `thermal_interface`, `k = 8.0`
  - `tungsten_insert`, `k = 170.0`
  - `air_gap`, `k = 0.026`

Later material regions in the PLC override earlier regions when a wedge centroid lies in more than one material box.

## Suggested large run

```bash
CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
mpirun --allow-run-as-root -np 4 ./build/heat_conduction_mpi \
  examples/large_composite_heat_multimaterial.plc i 110 90 8
```

The case intentionally omits the command-line conductivity argument so the PLC default conductivity and material regions are used.

## Benchmark record

The benchmarked `110 90 8` mesh has more than 100,000 wedge elements. Detailed process-count comparisons are recorded under `logs/benchmarks/` and summarized in the daily update log.

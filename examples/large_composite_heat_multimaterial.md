# Large composite multi-material heat-conduction case

Case file: `examples/large_composite_heat_multimaterial.plc`

This case is intended for stress-testing the heat-conduction pipeline above
100,000 wedge elements. The 2D footprint is a stitched, non-rectangular
polygonal silhouette that represents a composite made from several arbitrary
polygonal lobes. It is extruded through a thickness of `2.40` in `z`.

## Geometry and materials

- 22-vertex concave outer footprint
- extrusion height `H = 2.40`
- default matrix conductivity `k = 7.5`
- 12 centroid-assigned material regions:
  - `copper_bus_left`, `k = 385.0`
  - `copper_bus_right`, `k = 385.0`
  - `aluminum_frame`, `k = 205.0`
  - `ceramic_cap_upper`, `k = 24.0`
  - `graphite_midplane`, `k = 130.0`
  - `insulation_notch`, `k = 0.06`
  - `insulation_bridge`, `k = 0.08`
  - `steel_fastener_a`, `k = 16.0`
  - `steel_fastener_b`, `k = 16.0`
  - `diamond_spreader`, `k = 900.0`
  - `polymer_skin_top`, `k = 0.22`
  - `coolant_channel`, `k = 0.60`

Later material regions in the PLC override earlier ones for wedges whose
centroids lie in more than one material box.

## Benchmark command

```bash
CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

for np in 1 2 4; do
  /usr/bin/time -v mpirun --allow-run-as-root -np "${np}" \
    ./build/heat_conduction_mpi examples/large_composite_heat_multimaterial.plc i 110 90 10
done
```

The command intentionally omits the final conductivity override argument so the
PLC default and material-region conductivities are used.

## Verification and benchmark record

Actual results for the current branch are recorded in
`logs/benchmarks/large_composite_heat_multimaterial_2026-07-14.md` after the
build and MPI process-count sweep.

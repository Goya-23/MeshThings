# Complex multi-material heat-conduction case

Case file: `examples/complex_heat_multimaterial.plc`

This case exercises a larger non-rectangular 2D polygon with z extrusion and multiple material regions:

- 13-vertex concave footprint, extruded to height `1.60`
- default matrix conductivity `k = 4.0`
- rectangular material regions assigned by wedge centroid:
  - `copper_core`, `k = 385.0`
  - `ceramic_cap`, `k = 24.0`
  - `insulation_band`, `k = 0.08`
  - `graphite_midplane`, `k = 130.0`
  - `aluminum_lip`, `k = 205.0`

Later material regions in the PLC override earlier ones for wedges whose centroids lie in more than one material box.

## Suggested run

```bash
CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
mpirun --allow-run-as-root -np 2 ./build/heat_conduction_mpi \
  examples/complex_heat_multimaterial.plc i 18 18 6
```

The case intentionally omits the command-line conductivity argument so the PLC default conductivity and material regions are used.

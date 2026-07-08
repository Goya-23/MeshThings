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

## Verification record

Run on 2026-07-08:

```text
mpirun --allow-run-as-root -np 2 ./build/heat_conduction_mpi examples/complex_heat_multimaterial.plc i 18 18 6

Footprint CDT: 302 vertices, 451 triangles
Wedge mesh: 2114 vertices, 2706 wedge cells
OpenVolumeMesh cells: 2706, faces: 8059, vertices: 2114
Material default k = 4
  material[0] default: k = 4, wedge cells = 1147
  material[1] copper_core: k = 385, wedge cells = 378
  material[2] ceramic_cap: k = 24, wedge cells = 458
  material[3] insulation_band: k = 0.08, wedge cells = 258
  material[4] graphite_midplane: k = 130, wedge cells = 312
  material[5] aluminum_lip: k = 205, wedge cells = 153
METIS nodal partition into 2 parts
HYPRE PCG iterations: 24, residual: 2.7354e-11
HYPRE solve complete. Relative L2 error vs manufactured solution: 0.305397
```

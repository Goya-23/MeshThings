# Large composite heat-conduction benchmark - 2026-07-14

Case: `examples/large_composite_heat_multimaterial.plc`

Resolution: `nx = 110`, `ny = 90`, `nz = 10`

## Commands

```bash
CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

for np in 1 2 4; do
  /usr/bin/time -v mpirun --allow-run-as-root -np "${np}" \
    ./build/heat_conduction_mpi examples/large_composite_heat_multimaterial.plc i 110 90 10
done
```

## Results

Pending build and benchmark run on the current branch.

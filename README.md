# MeshThings

A mesh project demonstrating a full 2D heat conduction pipeline:

1. **CDT** — Constrained Delaunay triangulation of a PLC domain
2. **OpenVolumeMesh** — Surface mesh storage (`GeometricPolyhedralMeshV2d`)
3. **METIS** — Dynamic nodal mesh partitioning for MPI
4. **HYPRE** — Each rank assembles its sub-matrix rows; HYPRE PCG + BoomerAMG solves the gathered system on rank 0, then broadcasts the solution

## Problem

Steady-state heat conduction on the unit square:

\[
-\nabla \cdot (k \nabla u) = f
\]

with \(k = 1\), manufactured solution \(u = \sin(\pi x)\sin(\pi y)\), and Dirichlet boundary data from the exact solution.

## Dependencies

- CMake >= 3.16
- C++17 compiler
- MPI (OpenMPI)
- METIS (`libmetis-dev`)
- HYPRE (`libhypre-dev`)
- OpenVolumeMesh (fetched automatically via CMake `FetchContent`)

## Build

```bash
CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## Run

```bash
# 4 MPI ranks, incremental CDT, 16x16 interior grid
mpirun -np 4 ./build/heat_conduction_mpi examples/unit_square.plc i 16 16

# Use built-in unit square PLC
mpirun -np 2 ./build/heat_conduction_mpi
```

Arguments:

```
heat_conduction_mpi [plc_file] [method] [nx] [ny]
  plc_file  - optional PLC input (default: unit square)
  method    - i/I incremental CDT, s/S sweep-line alias (default: i)
  nx, ny    - interior Steiner grid resolution (default: 12)
```

PLC format (`examples/unit_square.plc`):

```
<number_of_boundary_vertices>
x y
...
<number_of_hole_seed_points>
x y
...
```

## Pipeline overview

```
PLC -> CDT -> OpenVolumeMesh -> METIS (rank 0) -> MPI broadcast
  -> per-rank FEM sub-matrix assembly -> HYPRE solve -> L2 error check
```

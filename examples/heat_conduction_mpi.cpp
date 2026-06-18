#include <cmath>
#include <iostream>
#include <memory>
#include <string>

#include <mpi.h>

#include <HYPRE.h>
#include <HYPRE_utilities.h>

#include "concreteDelaunayTriangulation.h"
#include "delaunayTriangulationFactory.h"
#include "heatFEM.h"
#include "hypreSolver.h"
#include "meshPartitioner.h"
#include "openVolumeMeshAdapter.h"
#include "plcParser.h"

namespace {

double computeL2Error(
    const Triangulation2D& mesh,
    const std::vector<double>& global_solution,
    int rank) {
    double local_error = 0.0;
    double local_area = 0.0;

    for (const auto& triangle : mesh.triangles) {
        const auto& p0 = mesh.vertices[triangle[0]];
        const auto& p1 = mesh.vertices[triangle[1]];
        const auto& p2 = mesh.vertices[triangle[2]];
        const double area2 = std::abs((p1.x_ - p0.x_) * (p2.y_ - p0.y_) - (p2.x_ - p0.x_) * (p1.y_ - p0.y_));
        const double area = 0.5 * area2;

        const double u0 = global_solution[triangle[0]];
        const double u1 = global_solution[triangle[1]];
        const double u2 = global_solution[triangle[2]];
        const double approx = (u0 + u1 + u2) / 3.0;

        const double x = (p0.x_ + p1.x_ + p2.x_) / 3.0;
        const double y = (p0.y_ + p1.y_ + p2.y_) / 3.0;
        const double exact = std::sin(M_PI * x) * std::sin(M_PI * y);
        const double diff = approx - exact;
        local_error += diff * diff * area;
        local_area += area;
    }

    double global_error = 0.0;
    double global_area = 0.0;
    MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_area, &global_area, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        return std::sqrt(global_error / global_area);
    }
    return 0.0;
}

}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int nx = 12;
    int ny = 12;
    char method = 'i';
    const char* plc_path = nullptr;

    if (argc >= 2) {
        plc_path = argv[1];
    }
    if (argc >= 3) {
        method = argv[2][0];
    }
    if (argc >= 5) {
        nx = std::stoi(argv[3]);
        ny = std::stoi(argv[4]);
    }

    try {
        HYPRE_Init();

        PLCParser parser;
        std::shared_ptr<PLC2D> plc = plc_path ? parser.parse(plc_path) : parser.makeUnitSquare();

        auto factory = std::make_shared<DelaunayTriangulationFactory>();
        std::shared_ptr<DelaunayTriangulation> triangulator = factory->produce(method);
        if (!triangulator) {
            throw std::runtime_error("Unknown triangulation method");
        }

        if (auto incremental = std::dynamic_pointer_cast<IncrementalDelaunayTriangulation>(triangulator)) {
            incremental->setGridDensity(nx, ny);
        } else if (auto sweep = std::dynamic_pointer_cast<SweepLineDelaunayTriangulation>(triangulator)) {
            sweep->setGridDensity(nx, ny);
        }

        triangulator->triangulate(plc);
        const Triangulation2D& triangulation = triangulator->result();

        SurfaceMesh2D surface_mesh = OpenVolumeMeshAdapter::buildSurfaceMesh(triangulation);
        Triangulation2D mesh_from_ovm = OpenVolumeMeshAdapter::extractTriangulation(surface_mesh);

        MeshPartition partition;
        if (rank == 0) {
            partition = MeshPartitioner::partitionNodal(mesh_from_ovm, size);
        }

        const int vertex_count = static_cast<int>(mesh_from_ovm.vertexCount());
        const int triangle_count = static_cast<int>(mesh_from_ovm.triangleCount());
        MPI_Bcast(&partition.num_parts, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (rank != 0) {
            partition.element_part.resize(triangle_count);
            partition.vertex_part.resize(vertex_count);
        }
        MPI_Bcast(partition.element_part.data(), triangle_count, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(partition.vertex_part.data(), vertex_count, MPI_INT, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            std::cout << "Mesh: " << mesh_from_ovm.vertexCount() << " vertices, "
                      << mesh_from_ovm.triangleCount() << " triangles\n";
            std::cout << "OpenVolumeMesh faces: " << surface_mesh.n_faces()
                      << ", vertices: " << surface_mesh.n_vertices() << "\n";
            std::cout << "METIS nodal partition into " << partition.num_parts << " parts\n";
        }

        mesh_from_ovm = MeshPartitioner::reorderByPartition(mesh_from_ovm, partition);

        LocalLinearSystem system = HeatFEMAssembler::assemble(MPI_COMM_WORLD, mesh_from_ovm, partition);
        std::vector<double> local_solution = HypreSolver::solve(MPI_COMM_WORLD, system);
        std::vector<double> gathered_solution = HypreSolver::gatherSolution(MPI_COMM_WORLD, system, local_solution);

        const double l2_error = computeL2Error(mesh_from_ovm, gathered_solution, rank);
        if (rank == 0) {
            std::cout << "HYPRE solve complete. Relative L2 error vs manufactured solution: " << l2_error << "\n";
        }

        HYPRE_Finalize();
    } catch (const std::exception& ex) {
        if (rank == 0) {
            std::cerr << "Error: " << ex.what() << "\n";
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Finalize();
    return 0;
}

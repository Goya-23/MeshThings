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
#include "wedgeMesh.h"

namespace {

double computeL2Error(
    const WedgeMesh3D& mesh,
    const std::vector<double>& global_solution,
    int rank) {
    double local_error = 0.0;
    double local_volume = 0.0;

    const std::array<std::array<int, 4>, 3> tetrahedra = {{
        {{0, 1, 2, 5}},
        {{0, 1, 5, 3}},
        {{1, 4, 5, 3}},
    }};

    for (const auto& wedge : mesh.wedges) {
        for (const auto& local_tet : tetrahedra) {
            std::array<int, 4> tet = {
                wedge[local_tet[0]],
                wedge[local_tet[1]],
                wedge[local_tet[2]],
                wedge[local_tet[3]],
            };

            const auto& p0 = mesh.vertices[tet[0]];
            const auto& p1 = mesh.vertices[tet[1]];
            const auto& p2 = mesh.vertices[tet[2]];
            const auto& p3 = mesh.vertices[tet[3]];

            const double v6 =
                (p1.x_ - p0.x_) * ((p2.y_ - p0.y_) * (p3.z_ - p0.z_) - (p2.z_ - p0.z_) * (p3.y_ - p0.y_)) -
                (p1.y_ - p0.y_) * ((p2.x_ - p0.x_) * (p3.z_ - p0.z_) - (p2.z_ - p0.z_) * (p3.x_ - p0.x_)) +
                (p1.z_ - p0.z_) * ((p2.x_ - p0.x_) * (p3.y_ - p0.y_) - (p2.y_ - p0.y_) * (p3.x_ - p0.x_));
            const double volume = std::abs(v6) / 6.0;

            const double u0 = global_solution[tet[0]];
            const double u1 = global_solution[tet[1]];
            const double u2 = global_solution[tet[2]];
            const double u3 = global_solution[tet[3]];
            const double approx = 0.25 * (u0 + u1 + u2 + u3);

            const double x = 0.25 * (p0.x_ + p1.x_ + p2.x_ + p3.x_);
            const double y = 0.25 * (p0.y_ + p1.y_ + p2.y_ + p3.y_);
            const double z = 0.25 * (p0.z_ + p1.z_ + p2.z_ + p3.z_);
            const double exact = std::sin(M_PI * x) * std::sin(M_PI * y) * std::sin(M_PI * z);
            const double diff = approx - exact;
            local_error += diff * diff * volume;
            local_volume += volume;
        }
    }

    double global_error = 0.0;
    double global_volume = 0.0;
    MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_volume, &global_volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        return std::sqrt(global_error / global_volume);
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

    int nx = 8;
    int ny = 8;
    int nz = 4;
    double height = 1.0;
    double conductivity = 1.0;
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
    if (argc >= 6) {
        nz = std::stoi(argv[5]);
    }
    if (argc >= 7) {
        height = std::stod(argv[6]);
    }
    if (argc >= 8) {
        conductivity = std::stod(argv[7]);
    }

    try {
        HYPRE_Init();

        PLCParser parser;
        std::shared_ptr<PLC2D> plc = plc_path ? parser.parse(plc_path) : parser.makeUnitSquare();
        if (argc < 7) {
            height = plc->getExtrusionHeight();
        }

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
        const Triangulation2D& footprint = triangulator->result();
        WedgeMesh3D wedge_mesh = WedgeExtruder::extrude(footprint, nz, height, conductivity);

        VolumeMesh3D volume_mesh = OpenVolumeMeshAdapter::buildWedgeVolumeMesh(wedge_mesh);
        WedgeMesh3D mesh_from_ovm =
            OpenVolumeMeshAdapter::extractWedgeMesh(volume_mesh, wedge_mesh.conductivity);

        MeshPartition partition;
        if (rank == 0) {
            partition = MeshPartitioner::partitionNodal(mesh_from_ovm, size);
        }

        const int vertex_count = static_cast<int>(mesh_from_ovm.vertexCount());
        const int wedge_count = static_cast<int>(mesh_from_ovm.wedgeCount());
        MPI_Bcast(&partition.num_parts, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (rank != 0) {
            partition.element_part.resize(wedge_count);
            partition.vertex_part.resize(vertex_count);
        }
        MPI_Bcast(partition.element_part.data(), wedge_count, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(partition.vertex_part.data(), vertex_count, MPI_INT, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            std::cout << "Footprint CDT: " << footprint.vertexCount() << " vertices, "
                      << footprint.triangleCount() << " triangles\n";
            std::cout << "Wedge mesh: " << mesh_from_ovm.vertexCount() << " vertices, "
                      << mesh_from_ovm.wedgeCount() << " wedge cells\n";
            std::cout << "OpenVolumeMesh cells: " << volume_mesh.n_cells()
                      << ", faces: " << volume_mesh.n_faces()
                      << ", vertices: " << volume_mesh.n_vertices() << "\n";
            std::cout << "Material conductivity k = " << mesh_from_ovm.conductivity << "\n";
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

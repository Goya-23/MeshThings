#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

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

struct MemorySampleKb {
    long long rss = 0;
    long long peak_rss = 0;
};

double beginTimedStage(MPI_Comm comm) {
    MPI_Barrier(comm);
    return MPI_Wtime();
}

double endTimedStage(MPI_Comm comm, double start_time) {
    const double local_seconds = MPI_Wtime() - start_time;
    double max_seconds = 0.0;
    MPI_Allreduce(&local_seconds, &max_seconds, 1, MPI_DOUBLE, MPI_MAX, comm);
    return max_seconds;
}

MemorySampleKb readProcessMemoryKb() {
    std::ifstream status("/proc/self/status");
    MemorySampleKb sample;
    std::string key;
    while (status >> key) {
        if (key == "VmRSS:") {
            status >> sample.rss;
        } else if (key == "VmHWM:") {
            status >> sample.peak_rss;
        }
        std::string rest_of_line;
        std::getline(status, rest_of_line);
    }
    return sample;
}

MemorySampleKb maxProcessMemoryKb(MPI_Comm comm) {
    const MemorySampleKb local = readProcessMemoryKb();
    MemorySampleKb global;
    MPI_Allreduce(&local.rss, &global.rss, 1, MPI_LONG_LONG, MPI_MAX, comm);
    MPI_Allreduce(&local.peak_rss, &global.peak_rss, 1, MPI_LONG_LONG, MPI_MAX, comm);
    return global;
}

double kibToMiB(long long kib) {
    return static_cast<double>(kib) / 1024.0;
}

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
    bool conductivity_overridden = false;
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
        conductivity_overridden = true;
    }

    try {
        HYPRE_Init();
        const double total_start = beginTimedStage(MPI_COMM_WORLD);

        double stage_start = beginTimedStage(MPI_COMM_WORLD);
        PLCParser parser;
        std::shared_ptr<PLC2D> plc = plc_path ? parser.parse(plc_path) : parser.makeUnitSquare();
        if (argc < 7) {
            height = plc->getExtrusionHeight();
        }
        if (!conductivity_overridden) {
            conductivity = plc->getDefaultConductivity();
        }
        const double parse_s = endTimedStage(MPI_COMM_WORLD, stage_start);

        stage_start = beginTimedStage(MPI_COMM_WORLD);
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
        const double triangulate_s = endTimedStage(MPI_COMM_WORLD, stage_start);

        stage_start = beginTimedStage(MPI_COMM_WORLD);
        const std::vector<MaterialRegion2D> material_regions =
            conductivity_overridden ? std::vector<MaterialRegion2D>{} : plc->getMaterialRegions();
        WedgeMesh3D wedge_mesh = WedgeExtruder::extrude(footprint, nz, height, conductivity, material_regions);
        const double extrude_s = endTimedStage(MPI_COMM_WORLD, stage_start);

        stage_start = beginTimedStage(MPI_COMM_WORLD);
        VolumeMesh3D volume_mesh = OpenVolumeMeshAdapter::buildWedgeVolumeMesh(wedge_mesh);
        WedgeMesh3D mesh_from_ovm = OpenVolumeMeshAdapter::extractWedgeMesh(
            volume_mesh,
            wedge_mesh.conductivity,
            wedge_mesh.wedge_conductivity,
            wedge_mesh.wedge_material_id,
            wedge_mesh.material_names,
            wedge_mesh.boundary_vertices);
        const double ovm_s = endTimedStage(MPI_COMM_WORLD, stage_start);

        stage_start = beginTimedStage(MPI_COMM_WORLD);
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
        const double partition_s = endTimedStage(MPI_COMM_WORLD, stage_start);

        if (rank == 0) {
            std::cout << "Footprint CDT: " << footprint.vertexCount() << " vertices, "
                      << footprint.triangleCount() << " triangles\n";
            std::cout << "Wedge mesh: " << mesh_from_ovm.vertexCount() << " vertices, "
                      << mesh_from_ovm.wedgeCount() << " wedge cells\n";
            std::cout << "OpenVolumeMesh cells: " << volume_mesh.n_cells()
                      << ", faces: " << volume_mesh.n_faces()
                      << ", vertices: " << volume_mesh.n_vertices() << "\n";
            std::map<int, std::size_t> material_cell_count;
            std::map<int, double> material_conductivity;
            for (std::size_t wedge_id = 0; wedge_id < mesh_from_ovm.wedgeCount(); ++wedge_id) {
                const int material_id =
                    wedge_id < mesh_from_ovm.wedge_material_id.size() ? mesh_from_ovm.wedge_material_id[wedge_id] : 0;
                material_cell_count[material_id]++;
                material_conductivity[material_id] = mesh_from_ovm.conductivityForWedge(wedge_id);
            }
            std::cout << "Material default k = " << mesh_from_ovm.conductivity << "\n";
            for (const auto& entry : material_cell_count) {
                const int material_id = entry.first;
                const std::string material_name =
                    material_id >= 0 && material_id < static_cast<int>(mesh_from_ovm.material_names.size())
                        ? mesh_from_ovm.material_names[material_id]
                        : std::string("material_") + std::to_string(material_id);
                std::cout << "  material[" << material_id << "] " << material_name
                          << ": k = " << material_conductivity[material_id]
                          << ", wedge cells = " << entry.second << "\n";
            }
            std::cout << "METIS nodal partition into " << partition.num_parts << " parts\n";
        }

        stage_start = beginTimedStage(MPI_COMM_WORLD);
        mesh_from_ovm = MeshPartitioner::reorderByPartition(mesh_from_ovm, partition);
        const double reorder_s = endTimedStage(MPI_COMM_WORLD, stage_start);

        stage_start = beginTimedStage(MPI_COMM_WORLD);
        LocalLinearSystem system = HeatFEMAssembler::assemble(MPI_COMM_WORLD, mesh_from_ovm, partition);
        const double assemble_s = endTimedStage(MPI_COMM_WORLD, stage_start);

        const long long local_rows = static_cast<long long>(system.owned_vertices.size());
        const long long local_nnz = static_cast<long long>(system.col_indices.size());
        long long total_rows = 0;
        long long min_rows = 0;
        long long max_rows = 0;
        long long total_nnz = 0;
        long long min_nnz = 0;
        long long max_nnz = 0;
        MPI_Reduce(&local_rows, &total_rows, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&local_rows, &min_rows, 1, MPI_LONG_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&local_rows, &max_rows, 1, MPI_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&local_nnz, &total_nnz, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&local_nnz, &min_nnz, 1, MPI_LONG_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&local_nnz, &max_nnz, 1, MPI_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            std::cout << "Sparse matrix rows per rank: min = " << min_rows
                      << ", max = " << max_rows
                      << ", total = " << total_rows << "\n";
            std::cout << "Sparse matrix nnz per rank: min = " << min_nnz
                      << ", max = " << max_nnz
                      << ", total = " << total_nnz << "\n";
        }

        stage_start = beginTimedStage(MPI_COMM_WORLD);
        std::vector<double> local_solution = HypreSolver::solve(MPI_COMM_WORLD, system);
        const double solve_s = endTimedStage(MPI_COMM_WORLD, stage_start);

        stage_start = beginTimedStage(MPI_COMM_WORLD);
        std::vector<double> gathered_solution = HypreSolver::gatherSolution(MPI_COMM_WORLD, system, local_solution);

        const double l2_error = computeL2Error(mesh_from_ovm, gathered_solution, rank);
        const double error_s = endTimedStage(MPI_COMM_WORLD, stage_start);
        const double total_s = endTimedStage(MPI_COMM_WORLD, total_start);
        const MemorySampleKb memory = maxProcessMemoryKb(MPI_COMM_WORLD);

        if (rank == 0) {
            std::cout << "HYPRE solve complete. Relative L2 error vs manufactured solution: " << l2_error << "\n";
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "Stage timings (max across ranks, seconds):"
                      << " parse=" << parse_s
                      << " triangulate=" << triangulate_s
                      << " extrude=" << extrude_s
                      << " ovm=" << ovm_s
                      << " partition=" << partition_s
                      << " reorder=" << reorder_s
                      << " assemble=" << assemble_s
                      << " solve=" << solve_s
                      << " error=" << error_s
                      << " total=" << total_s << "\n";
            std::cout << "Memory high-water mark (max rank): current_rss_mib="
                      << kibToMiB(memory.rss)
                      << " peak_rss_mib=" << kibToMiB(memory.peak_rss) << "\n";
            std::cout << "BENCHMARK"
                      << " case=" << (plc_path ? plc_path : "unit_square")
                      << " method=" << method
                      << " ranks=" << size
                      << " nx=" << nx
                      << " ny=" << ny
                      << " nz=" << nz
                      << " height=" << height
                      << " default_k=" << conductivity
                      << " footprint_vertices=" << footprint.vertexCount()
                      << " footprint_triangles=" << footprint.triangleCount()
                      << " wedge_vertices=" << mesh_from_ovm.vertexCount()
                      << " wedge_cells=" << mesh_from_ovm.wedgeCount()
                      << " matrix_rows=" << total_rows
                      << " matrix_nnz=" << total_nnz
                      << " l2_error=" << l2_error
                      << " parse_s=" << parse_s
                      << " triangulate_s=" << triangulate_s
                      << " extrude_s=" << extrude_s
                      << " ovm_s=" << ovm_s
                      << " partition_s=" << partition_s
                      << " reorder_s=" << reorder_s
                      << " assemble_s=" << assemble_s
                      << " solve_s=" << solve_s
                      << " error_s=" << error_s
                      << " total_s=" << total_s
                      << " rss_mib=" << kibToMiB(memory.rss)
                      << " peak_rss_mib=" << kibToMiB(memory.peak_rss)
                      << "\n";
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

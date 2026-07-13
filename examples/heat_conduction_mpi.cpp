#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
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

using Clock = std::chrono::steady_clock;

struct DoubleSummary {
    double min = 0.0;
    double avg = 0.0;
    double max = 0.0;
};

struct LongLongSummary {
    long long min = 0;
    double avg = 0.0;
    long long max = 0;
    long long sum = 0;
};

struct ProcessMemory {
    long long rss_kb = 0;
    long long peak_rss_kb = 0;
};

double secondsBetween(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double>(end - start).count();
}

DoubleSummary summarizeDouble(MPI_Comm comm, double local_value) {
    int size = 1;
    MPI_Comm_size(comm, &size);

    DoubleSummary summary;
    double sum = 0.0;
    MPI_Reduce(&local_value, &summary.min, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
    MPI_Reduce(&local_value, &summary.max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&local_value, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    summary.avg = sum / static_cast<double>(size);
    return summary;
}

LongLongSummary summarizeLongLong(MPI_Comm comm, long long local_value) {
    int size = 1;
    MPI_Comm_size(comm, &size);

    LongLongSummary summary;
    MPI_Reduce(&local_value, &summary.min, 1, MPI_LONG_LONG, MPI_MIN, 0, comm);
    MPI_Reduce(&local_value, &summary.max, 1, MPI_LONG_LONG, MPI_MAX, 0, comm);
    MPI_Reduce(&local_value, &summary.sum, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);
    summary.avg = static_cast<double>(summary.sum) / static_cast<double>(size);
    return summary;
}

void printStageTiming(MPI_Comm comm, int rank, const std::string& stage, double local_seconds) {
    const DoubleSummary timing = summarizeDouble(comm, local_seconds);
    if (rank == 0) {
        std::cout << "TIMING " << stage << " seconds min/avg/max = "
                  << timing.min << " / " << timing.avg << " / " << timing.max << "\n";
    }
}

ProcessMemory readProcessMemory() {
    ProcessMemory memory;
    std::ifstream status("/proc/self/status");
    std::string key;
    while (status >> key) {
        if (key == "VmRSS:" || key == "VmHWM:") {
            long long value = 0;
            std::string unit;
            status >> value >> unit;
            if (key == "VmRSS:") {
                memory.rss_kb = value;
            } else {
                memory.peak_rss_kb = value;
            }
        }
        std::string ignored;
        std::getline(status, ignored);
    }
    return memory;
}

void printMemorySummary(MPI_Comm comm, int rank, const ProcessMemory& local_memory) {
    const LongLongSummary rss = summarizeLongLong(comm, local_memory.rss_kb);
    const LongLongSummary peak = summarizeLongLong(comm, local_memory.peak_rss_kb);
    if (rank == 0) {
        std::cout << "MEMORY rss_mb avg/max = " << rss.avg / 1024.0 << " / " << rss.max / 1024.0
                  << ", peak_rss_mb avg/max = " << peak.avg / 1024.0 << " / " << peak.max / 1024.0
                  << "\n";
    }
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
        MPI_Barrier(MPI_COMM_WORLD);

        std::cout << std::setprecision(6);
        const auto total_start = Clock::now();
        auto stage_start = total_start;

        PLCParser parser;
        std::shared_ptr<PLC2D> plc = plc_path ? parser.parse(plc_path) : parser.makeUnitSquare();
        if (argc < 7) {
            height = plc->getExtrusionHeight();
        }
        if (!conductivity_overridden) {
            conductivity = plc->getDefaultConductivity();
        }
        printStageTiming(MPI_COMM_WORLD, rank, "parse_plc", secondsBetween(stage_start, Clock::now()));
        stage_start = Clock::now();

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
        printStageTiming(MPI_COMM_WORLD, rank, "triangulate_footprint", secondsBetween(stage_start, Clock::now()));
        stage_start = Clock::now();

        const std::vector<MaterialRegion2D> material_regions =
            conductivity_overridden ? std::vector<MaterialRegion2D>{} : plc->getMaterialRegions();
        WedgeMesh3D wedge_mesh = WedgeExtruder::extrude(footprint, nz, height, conductivity, material_regions);
        printStageTiming(MPI_COMM_WORLD, rank, "extrude_wedges", secondsBetween(stage_start, Clock::now()));
        stage_start = Clock::now();

        VolumeMesh3D volume_mesh = OpenVolumeMeshAdapter::buildWedgeVolumeMesh(wedge_mesh);
        WedgeMesh3D mesh_from_ovm = OpenVolumeMeshAdapter::extractWedgeMesh(
            volume_mesh,
            wedge_mesh.conductivity,
            wedge_mesh.wedge_conductivity,
            wedge_mesh.wedge_material_id,
            wedge_mesh.material_names,
            wedge_mesh.boundary_vertices);
        printStageTiming(MPI_COMM_WORLD, rank, "openvolumemesh_roundtrip", secondsBetween(stage_start, Clock::now()));
        stage_start = Clock::now();

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
        printStageTiming(MPI_COMM_WORLD, rank, "partition_and_broadcast", secondsBetween(stage_start, Clock::now()));
        stage_start = Clock::now();

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

        mesh_from_ovm = MeshPartitioner::reorderByPartition(mesh_from_ovm, partition);
        printStageTiming(MPI_COMM_WORLD, rank, "reorder_by_partition", secondsBetween(stage_start, Clock::now()));
        stage_start = Clock::now();

        LocalLinearSystem system = HeatFEMAssembler::assemble(MPI_COMM_WORLD, mesh_from_ovm, partition);
        printStageTiming(MPI_COMM_WORLD, rank, "assemble_local_system", secondsBetween(stage_start, Clock::now()));
        const LongLongSummary row_balance =
            summarizeLongLong(MPI_COMM_WORLD, static_cast<long long>(system.owned_vertices.size()));
        const LongLongSummary nnz_balance =
            summarizeLongLong(MPI_COMM_WORLD, static_cast<long long>(system.col_indices.size()));
        if (rank == 0) {
            std::cout << "MATRIX local_rows min/avg/max = "
                      << row_balance.min << " / " << row_balance.avg << " / " << row_balance.max
                      << ", local_nnz min/avg/max = "
                      << nnz_balance.min << " / " << nnz_balance.avg << " / " << nnz_balance.max << "\n";
        }
        stage_start = Clock::now();

        std::vector<double> local_solution = HypreSolver::solve(MPI_COMM_WORLD, system);
        printStageTiming(MPI_COMM_WORLD, rank, "hypre_solve", secondsBetween(stage_start, Clock::now()));
        stage_start = Clock::now();

        std::vector<double> gathered_solution = HypreSolver::gatherSolution(MPI_COMM_WORLD, system, local_solution);
        printStageTiming(MPI_COMM_WORLD, rank, "solution_gather", secondsBetween(stage_start, Clock::now()));
        stage_start = Clock::now();

        const double l2_error = computeL2Error(mesh_from_ovm, gathered_solution, rank);
        printStageTiming(MPI_COMM_WORLD, rank, "l2_error", secondsBetween(stage_start, Clock::now()));
        const double total_seconds = secondsBetween(total_start, Clock::now());
        const DoubleSummary total_summary = summarizeDouble(MPI_COMM_WORLD, total_seconds);
        const ProcessMemory local_memory = readProcessMemory();
        const LongLongSummary peak_memory_summary = summarizeLongLong(MPI_COMM_WORLD, local_memory.peak_rss_kb);
        printMemorySummary(MPI_COMM_WORLD, rank, local_memory);
        if (rank == 0) {
            std::cout << "HYPRE solve complete. Relative L2 error vs manufactured solution: " << l2_error << "\n";
            std::cout << "BENCHMARK"
                      << " processes=" << size
                      << " footprint_vertices=" << footprint.vertexCount()
                      << " footprint_triangles=" << footprint.triangleCount()
                      << " mesh_vertices=" << mesh_from_ovm.vertexCount()
                      << " wedge_cells=" << mesh_from_ovm.wedgeCount()
                      << " matrix_rows=" << system.global_size
                      << " matrix_nnz=" << nnz_balance.sum
                      << " total_seconds=" << total_summary.max
                      << " peak_rss_mb_max=" << peak_memory_summary.max / 1024.0
                      << " l2_error=" << l2_error << "\n";
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

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
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

long readProcStatusKb(const std::string& key) {
    std::ifstream status("/proc/self/status");
    std::string line;
    while (std::getline(status, line)) {
        if (line.rfind(key, 0) != 0) {
            continue;
        }

        std::istringstream stream(line.substr(key.size()));
        long value_kb = 0;
        stream >> value_kb;
        return value_kb;
    }
    return 0;
}

double reduceMaxDouble(MPI_Comm comm, double value, int root = 0) {
    double reduced = 0.0;
    MPI_Reduce(&value, &reduced, 1, MPI_DOUBLE, MPI_MAX, root, comm);
    return reduced;
}

long reduceMaxLong(MPI_Comm comm, long value, int root = 0) {
    long reduced = 0;
    MPI_Reduce(&value, &reduced, 1, MPI_LONG, MPI_MAX, root, comm);
    return reduced;
}

long long reduceSumLongLong(MPI_Comm comm, long long value, int root = 0) {
    long long reduced = 0;
    MPI_Reduce(&value, &reduced, 1, MPI_LONG_LONG, MPI_SUM, root, comm);
    return reduced;
}

int reduceMinInt(MPI_Comm comm, int value, int root = 0) {
    int reduced = std::numeric_limits<int>::max();
    MPI_Reduce(&value, &reduced, 1, MPI_INT, MPI_MIN, root, comm);
    return reduced;
}

int reduceMaxInt(MPI_Comm comm, int value, int root = 0) {
    int reduced = 0;
    MPI_Reduce(&value, &reduced, 1, MPI_INT, MPI_MAX, root, comm);
    return reduced;
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

        const double total_start = MPI_Wtime();

        double stage_start = MPI_Wtime();
        PLCParser parser;
        std::shared_ptr<PLC2D> plc = plc_path ? parser.parse(plc_path) : parser.makeUnitSquare();
        if (argc < 7) {
            height = plc->getExtrusionHeight();
        }
        if (!conductivity_overridden) {
            conductivity = plc->getDefaultConductivity();
        }
        const double parse_s = MPI_Wtime() - stage_start;

        stage_start = MPI_Wtime();
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
        const double cdt_s = MPI_Wtime() - stage_start;

        const std::vector<MaterialRegion2D> material_regions =
            conductivity_overridden ? std::vector<MaterialRegion2D>{} : plc->getMaterialRegions();

        stage_start = MPI_Wtime();
        WedgeMesh3D wedge_mesh = WedgeExtruder::extrude(footprint, nz, height, conductivity, material_regions);
        const double extrude_s = MPI_Wtime() - stage_start;

        stage_start = MPI_Wtime();
        VolumeMesh3D volume_mesh = OpenVolumeMeshAdapter::buildWedgeVolumeMesh(wedge_mesh);
        WedgeMesh3D mesh_from_ovm = OpenVolumeMeshAdapter::extractWedgeMesh(
            volume_mesh,
            wedge_mesh.conductivity,
            wedge_mesh.wedge_conductivity,
            wedge_mesh.wedge_material_id,
            wedge_mesh.material_names,
            wedge_mesh.boundary_vertices);
        const double ovm_s = MPI_Wtime() - stage_start;

        stage_start = MPI_Wtime();
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
        const double partition_s = MPI_Wtime() - stage_start;

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

        stage_start = MPI_Wtime();
        mesh_from_ovm = MeshPartitioner::reorderByPartition(mesh_from_ovm, partition);
        const double reorder_s = MPI_Wtime() - stage_start;

        stage_start = MPI_Wtime();
        LocalLinearSystem system = HeatFEMAssembler::assemble(MPI_COMM_WORLD, mesh_from_ovm, partition);
        const double assembly_s = MPI_Wtime() - stage_start;

        const int local_rows = static_cast<int>(system.owned_vertices.size());
        const int min_rows = reduceMinInt(MPI_COMM_WORLD, local_rows);
        const int max_rows = reduceMaxInt(MPI_COMM_WORLD, local_rows);
        const long long total_nnz = reduceSumLongLong(MPI_COMM_WORLD, static_cast<long long>(system.values.size()));
        const int local_nnz = static_cast<int>(system.values.size());
        const int min_nnz = reduceMinInt(MPI_COMM_WORLD, local_nnz);
        const int max_nnz = reduceMaxInt(MPI_COMM_WORLD, local_nnz);

        if (rank == 0) {
            std::cout << "Local rows per rank: min = " << min_rows
                      << ", max = " << max_rows
                      << ", total = " << system.global_size << "\n";
            std::cout << "Sparse matrix nnz: total = " << total_nnz
                      << ", min/rank = " << min_nnz
                      << ", max/rank = " << max_nnz << "\n";
        }

        stage_start = MPI_Wtime();
        std::vector<double> local_solution = HypreSolver::solve(MPI_COMM_WORLD, system);
        const double solve_s = MPI_Wtime() - stage_start;

        stage_start = MPI_Wtime();
        std::vector<double> gathered_solution = HypreSolver::gatherSolution(MPI_COMM_WORLD, system, local_solution);
        const double gather_solution_s = MPI_Wtime() - stage_start;

        stage_start = MPI_Wtime();
        const double l2_error = computeL2Error(mesh_from_ovm, gathered_solution, rank);
        const double error_s = MPI_Wtime() - stage_start;
        const double total_s = MPI_Wtime() - total_start;

        const double max_parse_s = reduceMaxDouble(MPI_COMM_WORLD, parse_s);
        const double max_cdt_s = reduceMaxDouble(MPI_COMM_WORLD, cdt_s);
        const double max_extrude_s = reduceMaxDouble(MPI_COMM_WORLD, extrude_s);
        const double max_ovm_s = reduceMaxDouble(MPI_COMM_WORLD, ovm_s);
        const double max_partition_s = reduceMaxDouble(MPI_COMM_WORLD, partition_s);
        const double max_reorder_s = reduceMaxDouble(MPI_COMM_WORLD, reorder_s);
        const double max_assembly_s = reduceMaxDouble(MPI_COMM_WORLD, assembly_s);
        const double max_solve_s = reduceMaxDouble(MPI_COMM_WORLD, solve_s);
        const double max_gather_solution_s = reduceMaxDouble(MPI_COMM_WORLD, gather_solution_s);
        const double max_error_s = reduceMaxDouble(MPI_COMM_WORLD, error_s);
        const double max_total_s = reduceMaxDouble(MPI_COMM_WORLD, total_s);
        const long max_rss_kb = reduceMaxLong(MPI_COMM_WORLD, readProcStatusKb("VmRSS:"));
        const long max_peak_kb = reduceMaxLong(MPI_COMM_WORLD, readProcStatusKb("VmHWM:"));

        if (rank == 0) {
            std::cout << "HYPRE solve complete. Relative L2 error vs manufactured solution: " << l2_error << "\n";
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "Stage timings (max seconds across ranks):"
                      << " parse=" << max_parse_s
                      << " cdt=" << max_cdt_s
                      << " extrude=" << max_extrude_s
                      << " ovm=" << max_ovm_s
                      << " partition=" << max_partition_s
                      << " reorder=" << max_reorder_s
                      << " assembly=" << max_assembly_s
                      << " solve=" << max_solve_s
                      << " gather_solution=" << max_gather_solution_s
                      << " l2=" << max_error_s
                      << " total=" << max_total_s << "\n";
            std::cout << "Memory usage (max across ranks): rss_mb="
                      << static_cast<double>(max_rss_kb) / 1024.0
                      << " peak_rss_mb=" << static_cast<double>(max_peak_kb) / 1024.0 << "\n";
            std::cout << "BENCHMARK"
                      << " ranks=" << size
                      << " plc=" << (plc_path ? plc_path : "unit_square")
                      << " method=" << method
                      << " nx=" << nx
                      << " ny=" << ny
                      << " nz=" << nz
                      << " height=" << height
                      << " footprint_vertices=" << footprint.vertexCount()
                      << " footprint_triangles=" << footprint.triangleCount()
                      << " vertices=" << mesh_from_ovm.vertexCount()
                      << " wedges=" << mesh_from_ovm.wedgeCount()
                      << " nnz=" << total_nnz
                      << " min_rows=" << min_rows
                      << " max_rows=" << max_rows
                      << " min_nnz=" << min_nnz
                      << " max_nnz=" << max_nnz
                      << " parse_s=" << max_parse_s
                      << " cdt_s=" << max_cdt_s
                      << " extrude_s=" << max_extrude_s
                      << " ovm_s=" << max_ovm_s
                      << " partition_s=" << max_partition_s
                      << " reorder_s=" << max_reorder_s
                      << " assembly_s=" << max_assembly_s
                      << " solve_s=" << max_solve_s
                      << " gather_solution_s=" << max_gather_solution_s
                      << " l2_s=" << max_error_s
                      << " total_s=" << max_total_s
                      << " rss_mb=" << static_cast<double>(max_rss_kb) / 1024.0
                      << " peak_rss_mb=" << static_cast<double>(max_peak_kb) / 1024.0
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

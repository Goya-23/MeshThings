#include "heatFEM.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <stdexcept>

namespace {

constexpr double kPi = 3.14159265358979323846;

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

Vec3 operator-(const Vec3& a, const Vec3& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    };
}

Vec3 operator/(const Vec3& vector, double scalar) {
    return {vector.x / scalar, vector.y / scalar, vector.z / scalar};
}

double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 toVec3(const Point3D& point) {
    return {point.x_, point.y_, point.z_};
}

}  // namespace

double HeatFEMAssembler::exactSolution(double x, double y, double z) {
    return std::sin(kPi * x) * std::sin(kPi * y) * std::sin(kPi * z);
}

double HeatFEMAssembler::sourceTerm(double x, double y, double z) {
    return 3.0 * kPi * kPi * std::sin(kPi * x) * std::sin(kPi * y) * std::sin(kPi * z);
}

bool HeatFEMAssembler::isBoundaryVertex(const WedgeMesh3D& mesh, int vertex_id) {
    if (vertex_id >= 0 && vertex_id < static_cast<int>(mesh.boundary_vertices.size())) {
        return mesh.boundary_vertices[vertex_id] != 0;
    }

    const Point3D& point = mesh.vertices[vertex_id];
    const double tol = 1.0e-8;

    double xmin = mesh.vertices.front().x_;
    double xmax = mesh.vertices.front().x_;
    double ymin = mesh.vertices.front().y_;
    double ymax = mesh.vertices.front().y_;
    double zmin = mesh.vertices.front().z_;
    double zmax = mesh.vertices.front().z_;
    for (const auto& vertex : mesh.vertices) {
        xmin = std::min(xmin, vertex.x_);
        xmax = std::max(xmax, vertex.x_);
        ymin = std::min(ymin, vertex.y_);
        ymax = std::max(ymax, vertex.y_);
        zmin = std::min(zmin, vertex.z_);
        zmax = std::max(zmax, vertex.z_);
    }

    const auto onBoundary = [&](double value, double min_value, double max_value) {
        return std::abs(value - min_value) < tol || std::abs(value - max_value) < tol;
    };

    return onBoundary(point.x_, xmin, xmax) || onBoundary(point.y_, ymin, ymax) ||
           onBoundary(point.z_, zmin, zmax);
}

void HeatFEMAssembler::addTetrahedronContribution(
    const WedgeMesh3D& mesh,
    const std::array<int, 4>& tet,
    double conductivity,
    std::map<int, std::map<int, double>>& matrix,
    std::map<int, double>& rhs) {
    const Vec3 p[4] = {
        toVec3(mesh.vertices[tet[0]]),
        toVec3(mesh.vertices[tet[1]]),
        toVec3(mesh.vertices[tet[2]]),
        toVec3(mesh.vertices[tet[3]]),
    };

    const double volume6 = dot(p[1] - p[0], cross(p[2] - p[0], p[3] - p[0]));
    if (std::abs(volume6) < 1.0e-14) {
        return;
    }
    const double volume = volume6 / 6.0;

    const Vec3 grad[4] = {
        cross(p[2] - p[1], p[3] - p[1]) / volume6,
        cross(p[3] - p[0], p[2] - p[0]) / volume6,
        cross(p[1] - p[0], p[3] - p[0]) / volume6,
        cross(p[1] - p[0], p[2] - p[0]) / volume6,
    };

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            matrix[tet[i]][tet[j]] += conductivity * volume * dot(grad[i], grad[j]);
        }
    }

    const Vec3 centroid = {
        0.25 * (p[0].x + p[1].x + p[2].x + p[3].x),
        0.25 * (p[0].y + p[1].y + p[2].y + p[3].y),
        0.25 * (p[0].z + p[1].z + p[2].z + p[3].z),
    };
    const double load = conductivity * sourceTerm(centroid.x, centroid.y, centroid.z) * volume / 4.0;
    for (int i = 0; i < 4; ++i) {
        rhs[tet[i]] += load;
    }
}

void HeatFEMAssembler::addWedgeContribution(
    const WedgeMesh3D& mesh,
    std::size_t wedge_id,
    const std::array<int, 6>& wedge,
    std::map<int, std::map<int, double>>& matrix,
    std::map<int, double>& rhs) {
    const std::array<std::array<int, 4>, 3> tetrahedra = {{
        {{wedge[0], wedge[1], wedge[2], wedge[5]}},
        {{wedge[0], wedge[1], wedge[5], wedge[3]}},
        {{wedge[1], wedge[4], wedge[5], wedge[3]}},
    }};

    const double conductivity = mesh.conductivityForWedge(wedge_id);
    for (const auto& tet : tetrahedra) {
        addTetrahedronContribution(mesh, tet, conductivity, matrix, rhs);
    }
}

LocalLinearSystem HeatFEMAssembler::assemble(
    MPI_Comm comm,
    const WedgeMesh3D& mesh,
    const MeshPartition& partition) {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (static_cast<int>(partition.vertex_part.size()) != static_cast<int>(mesh.vertexCount())) {
        throw std::runtime_error("Partition size does not match mesh vertex count");
    }
    if (size != partition.num_parts) {
        throw std::runtime_error("MPI communicator size must match METIS partition count");
    }

    LocalLinearSystem system;
    system.global_size = static_cast<int>(mesh.vertexCount());

    std::vector<int> row_counts(size, 0);
    for (int vertex = 0; vertex < system.global_size; ++vertex) {
        if (partition.vertex_part[vertex] == rank) {
            ++row_counts[rank];
        }
    }

    std::vector<int> global_row_counts(size, 0);
    MPI_Allgather(&row_counts[rank], 1, MPI_INT, global_row_counts.data(), 1, MPI_INT, comm);

    std::vector<int> row_displs(size + 1, 0);
    for (int process = 0; process < size; ++process) {
        row_displs[process + 1] = row_displs[process] + global_row_counts[process];
    }

    system.row_start = row_displs[rank];
    system.row_end = row_displs[rank + 1] - 1;

    for (int vertex = 0; vertex < system.global_size; ++vertex) {
        if (partition.vertex_part[vertex] == rank) {
            system.owned_vertices.push_back(vertex);
        }
    }

    std::set<int> touched_vertices(system.owned_vertices.begin(), system.owned_vertices.end());
    std::vector<std::size_t> incident_wedges;
    incident_wedges.reserve(mesh.wedges.size() / static_cast<std::size_t>(size) + 1);
    for (std::size_t wedge_id = 0; wedge_id < mesh.wedges.size(); ++wedge_id) {
        const auto& wedge = mesh.wedges[wedge_id];
        bool touches_owned = false;
        for (int corner = 0; corner < 6; ++corner) {
            if (partition.vertex_part[wedge[corner]] == rank) {
                touches_owned = true;
                break;
            }
        }
        if (!touches_owned) {
            continue;
        }
        incident_wedges.push_back(wedge_id);
        for (int corner = 0; corner < 6; ++corner) {
            touched_vertices.insert(wedge[corner]);
        }
    }

    system.local_to_global.assign(touched_vertices.begin(), touched_vertices.end());
    for (std::size_t local = 0; local < system.local_to_global.size(); ++local) {
        system.global_to_local[system.local_to_global[local]] = static_cast<int>(local);
    }

    std::map<int, std::map<int, double>> global_matrix;
    std::map<int, double> global_rhs;
    for (std::size_t wedge_id : incident_wedges) {
        addWedgeContribution(mesh, wedge_id, mesh.wedges[wedge_id], global_matrix, global_rhs);
    }

    system.is_dirichlet.assign(system.local_to_global.size(), 0);
    for (std::size_t local = 0; local < system.local_to_global.size(); ++local) {
        const int global = system.local_to_global[local];
        if (isBoundaryVertex(mesh, global)) {
            system.is_dirichlet[local] = 1;
        }
    }

    const int local_rows = static_cast<int>(system.owned_vertices.size());
    system.row_ptr.assign(local_rows + 1, 0);
    system.row_indices.resize(local_rows);
    for (int local_row = 0; local_row < local_rows; ++local_row) {
        system.row_indices[local_row] = system.owned_vertices[local_row];
    }

    std::vector<std::vector<std::pair<int, double>>> rows(local_rows);
    system.rhs.assign(local_rows, 0.0);
    system.solution.assign(local_rows, 0.0);

    for (int local_row = 0; local_row < local_rows; ++local_row) {
        const int global_row = system.owned_vertices[local_row];
        const Point3D& point = mesh.vertices[global_row];
        if (system.is_dirichlet[system.global_to_local[global_row]]) {
            rows[local_row].push_back({global_row, 1.0});
            system.rhs[local_row] = exactSolution(point.x_, point.y_, point.z_);
            continue;
        }

        for (const auto& column_entry : global_matrix[global_row]) {
            const int global_col = column_entry.first;
            const double value = column_entry.second;
            const Point3D& col_point = mesh.vertices[global_col];
            if (isBoundaryVertex(mesh, global_col)) {
                system.rhs[local_row] -= value * exactSolution(col_point.x_, col_point.y_, col_point.z_);
            } else {
                rows[local_row].push_back({global_col, value});
            }
        }
        system.rhs[local_row] += global_rhs[global_row];
    }

    for (int local_row = 0; local_row < local_rows; ++local_row) {
        std::sort(rows[local_row].begin(), rows[local_row].end(), [](const auto& a, const auto& b) {
            return a.first < b.first;
        });
        system.row_ptr[local_row + 1] = system.row_ptr[local_row] + static_cast<int>(rows[local_row].size());
    }

    system.col_indices.clear();
    system.values.clear();
    system.col_indices.reserve(system.row_ptr.back());
    system.values.reserve(system.row_ptr.back());
    for (int local_row = 0; local_row < local_rows; ++local_row) {
        for (const auto& entry : rows[local_row]) {
            system.col_indices.push_back(entry.first);
            system.values.push_back(entry.second);
        }
    }

    return system;
}

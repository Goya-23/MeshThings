#include "heatFEM.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

namespace {

constexpr double kPi = 3.14159265358979323846;

bool isOnSegment(const Point2D& p, const Point2D& a, const Point2D& b) {
    const double cross = (b.x_ - a.x_) * (p.y_ - a.y_) - (b.y_ - a.y_) * (p.x_ - a.x_);
    if (std::abs(cross) > 1.0e-10) {
        return false;
    }
    const double dot = (p.x_ - a.x_) * (p.x_ - b.x_) + (p.y_ - a.y_) * (p.y_ - b.y_);
    return dot <= 1.0e-10;
}

}  // namespace

double HeatFEMAssembler::exactSolution(double x, double y) {
    return std::sin(kPi * x) * std::sin(kPi * y);
}

double HeatFEMAssembler::sourceTerm(double x, double y) {
    return 2.0 * kPi * kPi * std::sin(kPi * x) * std::sin(kPi * y);
}

bool HeatFEMAssembler::isBoundaryVertex(const Triangulation2D& mesh, int vertex_id) {
    const Point2D& point = mesh.vertices[vertex_id];
    const double tol = 1.0e-8;
    const auto onBoundary = [&](double value, double min_value, double max_value) {
        return std::abs(value - min_value) < tol || std::abs(value - max_value) < tol;
    };

    double xmin = mesh.vertices.front().x_;
    double xmax = mesh.vertices.front().x_;
    double ymin = mesh.vertices.front().y_;
    double ymax = mesh.vertices.front().y_;
    for (const auto& vertex : mesh.vertices) {
        xmin = std::min(xmin, vertex.x_);
        xmax = std::max(xmax, vertex.x_);
        ymin = std::min(ymin, vertex.y_);
        ymax = std::max(ymax, vertex.y_);
    }

    return onBoundary(point.x_, xmin, xmax) || onBoundary(point.y_, ymin, ymax);
}

void HeatFEMAssembler::addTriangleContribution(
    const Triangulation2D& mesh,
    const std::array<int, 3>& tri,
    double conductivity,
    std::map<int, std::map<int, double>>& matrix,
    std::map<int, double>& rhs) {
    const Point2D& p0 = mesh.vertices[tri[0]];
    const Point2D& p1 = mesh.vertices[tri[1]];
    const Point2D& p2 = mesh.vertices[tri[2]];

    const double area2 = std::abs((p1.x_ - p0.x_) * (p2.y_ - p0.y_) - (p2.x_ - p0.x_) * (p1.y_ - p0.y_));
    if (area2 < 1.0e-14) {
        return;
    }
    const double area = 0.5 * area2;

    const double b0 = p1.y_ - p2.y_;
    const double b1 = p2.y_ - p0.y_;
    const double b2 = p0.y_ - p1.y_;
    const double c0 = p2.x_ - p1.x_;
    const double c1 = p0.x_ - p2.x_;
    const double c2 = p1.x_ - p0.x_;

    const double factor = conductivity / area;
    const int ids[3] = {tri[0], tri[1], tri[2]};
    const double b[3] = {b0, b1, b2};
    const double c[3] = {c0, c1, c2};

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            matrix[ids[i]][ids[j]] += factor * (b[i] * b[j] + c[i] * c[j]);
        }
    }

    const double centroid_x = (p0.x_ + p1.x_ + p2.x_) / 3.0;
    const double centroid_y = (p0.y_ + p1.y_ + p2.y_) / 3.0;
    const double load = sourceTerm(centroid_x, centroid_y) * area / 3.0;
    for (int i = 0; i < 3; ++i) {
        rhs[ids[i]] += load;
    }
}

LocalLinearSystem HeatFEMAssembler::assemble(
    MPI_Comm comm,
    const Triangulation2D& mesh,
    const MeshPartition& partition,
    double conductivity) {
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
    for (std::size_t element = 0; element < mesh.triangleCount(); ++element) {
        const auto& tri = mesh.triangles[element];
        bool touches_owned = false;
        for (int corner = 0; corner < 3; ++corner) {
            if (partition.vertex_part[tri[corner]] == rank) {
                touches_owned = true;
                break;
            }
        }
        if (!touches_owned) {
            continue;
        }
        for (int corner = 0; corner < 3; ++corner) {
            touched_vertices.insert(tri[corner]);
        }
    }

    system.local_to_global.assign(touched_vertices.begin(), touched_vertices.end());
    for (std::size_t local = 0; local < system.local_to_global.size(); ++local) {
        system.global_to_local[system.local_to_global[local]] = static_cast<int>(local);
    }

    std::map<int, std::map<int, double>> global_matrix;
    std::map<int, double> global_rhs;
    for (const auto& triangle : mesh.triangles) {
        addTriangleContribution(mesh, triangle, conductivity, global_matrix, global_rhs);
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
        if (system.is_dirichlet[system.global_to_local[global_row]]) {
            rows[local_row].push_back({global_row, 1.0});
            system.rhs[local_row] = exactSolution(mesh.vertices[global_row].x_, mesh.vertices[global_row].y_);
            continue;
        }

        for (const auto& column_entry : global_matrix[global_row]) {
            const int global_col = column_entry.first;
            const double value = column_entry.second;
            if (isBoundaryVertex(mesh, global_col)) {
                system.rhs[local_row] -= value * exactSolution(
                    mesh.vertices[global_col].x_, mesh.vertices[global_col].y_);
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

#include "constrainedDelaunay.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

constexpr double kEpsilon = 1.0e-12;
constexpr int kStructuredClipCellThreshold = 4096;

int findVertexIndex(const Triangulation2D& mesh, const Point2D& point) {
    for (std::size_t i = 0; i < mesh.vertices.size(); ++i) {
        const double dx = mesh.vertices[i].x_ - point.x_;
        const double dy = mesh.vertices[i].y_ - point.y_;
        if ((dx * dx + dy * dy) < kEpsilon) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

int addVertex(Triangulation2D& mesh, const Point2D& point) {
    const int existing = findVertexIndex(mesh, point);
    if (existing >= 0) {
        return existing;
    }
    mesh.vertices.push_back(point);
    return static_cast<int>(mesh.vertices.size()) - 1;
}

double cross2d(const Point2D& a, const Point2D& b, const Point2D& c) {
    return (b.x_ - a.x_) * (c.y_ - a.y_) - (b.y_ - a.y_) * (c.x_ - a.x_);
}

bool pointInTriangle(const Point2D& p, const Point2D& a, const Point2D& b, const Point2D& c) {
    const double c0 = cross2d(a, b, p);
    const double c1 = cross2d(b, c, p);
    const double c2 = cross2d(c, a, p);
    const bool has_neg = (c0 < -kEpsilon) || (c1 < -kEpsilon) || (c2 < -kEpsilon);
    const bool has_pos = (c0 > kEpsilon) || (c1 > kEpsilon) || (c2 > kEpsilon);
    return !(has_neg && has_pos);
}

bool pointOnSegment(const Point2D& point, const Point2D& a, const Point2D& b) {
    if (std::abs(cross2d(a, b, point)) > kEpsilon) {
        return false;
    }
    return point.x_ >= std::min(a.x_, b.x_) - kEpsilon &&
           point.x_ <= std::max(a.x_, b.x_) + kEpsilon &&
           point.y_ >= std::min(a.y_, b.y_) - kEpsilon &&
           point.y_ <= std::max(a.y_, b.y_) + kEpsilon;
}

bool pointInPolygon(const Point2D& point, const std::vector<Point2D>& polygon) {
    bool inside = false;
    for (std::size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++) {
        const Point2D& a = polygon[i];
        const Point2D& b = polygon[j];
        if (pointOnSegment(point, a, b)) {
            return true;
        }
        const bool edge_crosses_ray = (a.y_ > point.y_) != (b.y_ > point.y_);
        if (edge_crosses_ray) {
            const double x_intersection =
                (b.x_ - a.x_) * (point.y_ - a.y_) / (b.y_ - a.y_) + a.x_;
            if (point.x_ < x_intersection) {
                inside = !inside;
            }
        }
    }
    return inside;
}

bool isAxisAlignedRectangle(const PLC2D& plc) {
    if (plc.getPolygons().empty()) {
        return false;
    }
    const auto& boundary = plc.getPolygons().front().getPoints();
    if (boundary.size() != 4) {
        return false;
    }

    int axis_aligned_edges = 0;
    for (std::size_t i = 0; i < boundary.size(); ++i) {
        const auto& a = boundary[i];
        const auto& b = boundary[(i + 1) % boundary.size()];
        if (std::abs(a.x_ - b.x_) < kEpsilon || std::abs(a.y_ - b.y_) < kEpsilon) {
            ++axis_aligned_edges;
        }
    }
    return axis_aligned_edges == 4;
}

std::pair<int, int> sortedEdge(int a, int b) {
    if (a > b) {
        std::swap(a, b);
    }
    return {a, b};
}

void compactUnusedVertices(Triangulation2D& mesh) {
    std::vector<char> used(mesh.vertices.size(), 0);
    for (const auto& triangle : mesh.triangles) {
        used[triangle[0]] = 1;
        used[triangle[1]] = 1;
        used[triangle[2]] = 1;
    }
    for (const auto& edge : mesh.constrained_edges) {
        if (edge.first >= 0 && edge.first < static_cast<int>(used.size())) {
            used[edge.first] = 1;
        }
        if (edge.second >= 0 && edge.second < static_cast<int>(used.size())) {
            used[edge.second] = 1;
        }
    }

    std::vector<int> old_to_new(mesh.vertices.size(), -1);
    std::vector<Point2D> vertices;
    vertices.reserve(mesh.vertices.size());
    for (std::size_t vertex_id = 0; vertex_id < mesh.vertices.size(); ++vertex_id) {
        if (!used[vertex_id]) {
            continue;
        }
        old_to_new[vertex_id] = static_cast<int>(vertices.size());
        vertices.push_back(mesh.vertices[vertex_id]);
    }

    for (auto& triangle : mesh.triangles) {
        triangle[0] = old_to_new[triangle[0]];
        triangle[1] = old_to_new[triangle[1]];
        triangle[2] = old_to_new[triangle[2]];
    }

    std::vector<std::pair<int, int>> constrained_edges;
    constrained_edges.reserve(mesh.constrained_edges.size());
    for (const auto& edge : mesh.constrained_edges) {
        if (edge.first < 0 || edge.second < 0 ||
            edge.first >= static_cast<int>(old_to_new.size()) ||
            edge.second >= static_cast<int>(old_to_new.size())) {
            continue;
        }
        const int first = old_to_new[edge.first];
        const int second = old_to_new[edge.second];
        if (first >= 0 && second >= 0 && first != second) {
            constrained_edges.emplace_back(first, second);
        }
    }

    mesh.vertices = std::move(vertices);
    mesh.constrained_edges = std::move(constrained_edges);
}

void triangulateStructuredRectangle(int nx, int ny, Triangulation2D& mesh) {
    if (mesh.vertices.empty()) {
        return;
    }

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

    nx = std::max(nx, 1);
    ny = std::max(ny, 1);
    mesh.vertices.clear();
    mesh.triangles.clear();

    const int npx = nx + 1;
    const int npy = ny + 1;
    mesh.vertices.reserve(static_cast<std::size_t>(npx * npy));

    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            const double x = xmin + (xmax - xmin) * static_cast<double>(i) / static_cast<double>(nx);
            const double y = ymin + (ymax - ymin) * static_cast<double>(j) / static_cast<double>(ny);
            mesh.vertices.emplace_back(x, y);
        }
    }

    auto vertex_index = [npx](int i, int j) { return j * npx + i; };
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const int v00 = vertex_index(i, j);
            const int v10 = vertex_index(i + 1, j);
            const int v01 = vertex_index(i, j + 1);
            const int v11 = vertex_index(i + 1, j + 1);
            mesh.triangles.push_back({v00, v10, v11});
            mesh.triangles.push_back({v00, v11, v01});
        }
    }

    mesh.constrained_edges.clear();
    for (int i = 0; i < nx; ++i) {
        mesh.constrained_edges.emplace_back(vertex_index(i, 0), vertex_index(i + 1, 0));
        mesh.constrained_edges.emplace_back(vertex_index(i, ny), vertex_index(i + 1, ny));
    }
    for (int j = 0; j < ny; ++j) {
        mesh.constrained_edges.emplace_back(vertex_index(0, j), vertex_index(0, j + 1));
        mesh.constrained_edges.emplace_back(vertex_index(nx, j), vertex_index(nx, j + 1));
    }
}

bool triangulateStructuredClip(const PLC2D& plc, int nx, int ny, Triangulation2D& mesh) {
    if (plc.getPolygons().empty()) {
        return false;
    }

    const auto& boundary = plc.getPolygons().front().getPoints();
    if (boundary.size() < 3) {
        return false;
    }

    double xmin = boundary.front().x_;
    double xmax = boundary.front().x_;
    double ymin = boundary.front().y_;
    double ymax = boundary.front().y_;
    for (const auto& point : boundary) {
        xmin = std::min(xmin, point.x_);
        xmax = std::max(xmax, point.x_);
        ymin = std::min(ymin, point.y_);
        ymax = std::max(ymax, point.y_);
    }

    nx = std::max(nx, 1);
    ny = std::max(ny, 1);
    const int npx = nx + 1;
    const int npy = ny + 1;
    mesh.vertices.clear();
    mesh.triangles.clear();
    mesh.constrained_edges.clear();
    mesh.vertices.reserve(static_cast<std::size_t>(npx) * static_cast<std::size_t>(npy));

    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            const double x = xmin + (xmax - xmin) * static_cast<double>(i) / static_cast<double>(nx);
            const double y = ymin + (ymax - ymin) * static_cast<double>(j) / static_cast<double>(ny);
            mesh.vertices.emplace_back(x, y);
        }
    }

    auto vertex_index = [npx](int i, int j) { return j * npx + i; };
    const auto keep_triangle = [&](const std::array<int, 3>& triangle) {
        const auto& a = mesh.vertices[triangle[0]];
        const auto& b = mesh.vertices[triangle[1]];
        const auto& c = mesh.vertices[triangle[2]];
        const Point2D centroid(
            (a.x_ + b.x_ + c.x_) / 3.0,
            (a.y_ + b.y_ + c.y_) / 3.0);
        if (!pointInPolygon(centroid, boundary)) {
            return false;
        }
        for (const auto& hole : plc.getHoles()) {
            if (pointInTriangle(hole, a, b, c)) {
                return false;
            }
        }
        return true;
    };

    mesh.triangles.reserve(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) * 2);
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const int v00 = vertex_index(i, j);
            const int v10 = vertex_index(i + 1, j);
            const int v01 = vertex_index(i, j + 1);
            const int v11 = vertex_index(i + 1, j + 1);
            const std::array<int, 3> lower = {v00, v10, v11};
            const std::array<int, 3> upper = {v00, v11, v01};
            if (keep_triangle(lower)) {
                mesh.triangles.push_back(lower);
            }
            if (keep_triangle(upper)) {
                mesh.triangles.push_back(upper);
            }
        }
    }

    if (mesh.triangles.empty()) {
        return false;
    }

    struct EdgeHash {
        std::size_t operator()(const std::pair<int, int>& edge) const {
            return std::hash<int>()(edge.first) ^ (std::hash<int>()(edge.second) << 1);
        }
    };
    std::unordered_map<std::pair<int, int>, int, EdgeHash> edge_count;
    for (const auto& triangle : mesh.triangles) {
        ++edge_count[sortedEdge(triangle[0], triangle[1])];
        ++edge_count[sortedEdge(triangle[1], triangle[2])];
        ++edge_count[sortedEdge(triangle[2], triangle[0])];
    }
    mesh.constrained_edges.reserve(edge_count.size());
    for (const auto& entry : edge_count) {
        if (entry.second == 1) {
            mesh.constrained_edges.push_back(entry.first);
        }
    }

    compactUnusedVertices(mesh);
    return true;
}

}  // namespace

Triangulation2D ConstrainedDelaunayTriangulator::triangulate(const PLC2D& plc, int nx, int ny) {
    Triangulation2D mesh;
    addInteriorPoints(plc, nx, ny, mesh);
    buildBoundaryConstraints(plc, mesh);

    if (isAxisAlignedRectangle(plc)) {
        triangulateStructuredRectangle(nx, ny, mesh);
        return mesh;
    }

    if (nx * ny >= kStructuredClipCellThreshold &&
        triangulateStructuredClip(plc, nx, ny, mesh)) {
        return mesh;
    }

    bowyerWatson(mesh);
    enforceConstraints(mesh);
    removeExteriorTriangles(plc, mesh);
    return mesh;
}

void ConstrainedDelaunayTriangulator::addInteriorPoints(
    const PLC2D& plc, int nx, int ny, Triangulation2D& mesh) {
    if (plc.getPolygons().empty()) {
        return;
    }

    const auto& boundary = plc.getPolygons().front().getPoints();
    for (const auto& point : boundary) {
        addVertex(mesh, point);
    }

    if (isAxisAlignedRectangle(plc)) {
        return;
    }

    double xmin = boundary.front().x_;
    double xmax = boundary.front().x_;
    double ymin = boundary.front().y_;
    double ymax = boundary.front().y_;
    for (const auto& point : boundary) {
        xmin = std::min(xmin, point.x_);
        xmax = std::max(xmax, point.x_);
        ymin = std::min(ymin, point.y_);
        ymax = std::max(ymax, point.y_);
    }

    nx = std::max(nx, 1);
    ny = std::max(ny, 1);
    for (int j = 1; j < ny; ++j) {
        for (int i = 1; i < nx; ++i) {
            const double x = xmin + (xmax - xmin) * static_cast<double>(i) / static_cast<double>(nx);
            const double y = ymin + (ymax - ymin) * static_cast<double>(j) / static_cast<double>(ny);
            addVertex(mesh, Point2D(x, y));
        }
    }
}

void ConstrainedDelaunayTriangulator::buildBoundaryConstraints(
    const PLC2D& plc, Triangulation2D& mesh) {
    mesh.constrained_edges.clear();
    if (plc.getPolygons().empty()) {
        return;
    }

    const auto& boundary = plc.getPolygons().front().getPoints();
    const int n = static_cast<int>(boundary.size());
    for (int i = 0; i < n; ++i) {
        const int v0 = findVertexIndex(mesh, boundary[i]);
        const int v1 = findVertexIndex(mesh, boundary[(i + 1) % n]);
        if (v0 >= 0 && v1 >= 0) {
            mesh.constrained_edges.emplace_back(v0, v1);
        }
    }
}

bool ConstrainedDelaunayTriangulator::inCircumcircle(
    const Point2D& p, const Point2D& a, const Point2D& b, const Point2D& c) {
    const double ax = a.x_ - p.x_;
    const double ay = a.y_ - p.y_;
    const double bx = b.x_ - p.x_;
    const double by = b.y_ - p.y_;
    const double cx = c.x_ - p.x_;
    const double cy = c.y_ - p.y_;
    const double det =
        (ax * ax + ay * ay) * (bx * cy - by * cx) -
        (bx * bx + by * by) * (ax * cy - ay * cx) +
        (cx * cx + cy * cy) * (ax * by - ay * bx);
    return det > kEpsilon;
}

void ConstrainedDelaunayTriangulator::bowyerWatson(Triangulation2D& mesh) {
    if (mesh.vertices.size() < 3) {
        return;
    }

    double xmin = mesh.vertices[0].x_;
    double xmax = mesh.vertices[0].x_;
    double ymin = mesh.vertices[0].y_;
    double ymax = mesh.vertices[0].y_;
    for (const auto& vertex : mesh.vertices) {
        xmin = std::min(xmin, vertex.x_);
        xmax = std::max(xmax, vertex.x_);
        ymin = std::min(ymin, vertex.y_);
        ymax = std::max(ymax, vertex.y_);
    }

    const double dx = xmax - xmin;
    const double dy = ymax - ymin;
    const double delta = std::max(dx, dy);
    const double scale = delta > 0.0 ? 20.0 * delta : 1.0;

    const int i0 = static_cast<int>(mesh.vertices.size());
    mesh.vertices.emplace_back(xmin - scale, ymin - scale);
    const int i1 = static_cast<int>(mesh.vertices.size());
    mesh.vertices.emplace_back(xmax + scale, ymin - scale);
    const int i2 = static_cast<int>(mesh.vertices.size());
    mesh.vertices.emplace_back(xmin + 0.5 * dx, ymax + scale);

    mesh.triangles.clear();
    mesh.triangles.push_back({i0, i1, i2});

    const int original_vertex_count = i0;
    for (int vi = 0; vi < original_vertex_count; ++vi) {
        const Point2D& point = mesh.vertices[vi];
        std::vector<int> bad_triangles;
        for (std::size_t ti = 0; ti < mesh.triangles.size(); ++ti) {
            const auto& tri = mesh.triangles[ti];
            if (inCircumcircle(point, mesh.vertices[tri[0]], mesh.vertices[tri[1]], mesh.vertices[tri[2]])) {
                bad_triangles.push_back(static_cast<int>(ti));
            }
        }

        struct EdgeHash {
            std::size_t operator()(const std::pair<int, int>& edge) const {
                return std::hash<int>()(edge.first) ^ (std::hash<int>()(edge.second) << 1);
            }
        };
        std::unordered_map<std::pair<int, int>, int, EdgeHash> edge_count;

        for (int bad_index : bad_triangles) {
            const auto& tri = mesh.triangles[bad_index];
            for (int k = 0; k < 3; ++k) {
                const int a = tri[k];
                const int b = tri[(k + 1) % 3];
                const auto edge = (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
                edge_count[edge]++;
            }
        }

        std::vector<std::pair<int, int>> polygon;
        for (const auto& entry : edge_count) {
            if (entry.second == 1) {
                polygon.push_back(entry.first);
            }
        }

        std::vector<std::array<int, 3>> remaining;
        std::unordered_set<int> bad_set(bad_triangles.begin(), bad_triangles.end());
        for (std::size_t ti = 0; ti < mesh.triangles.size(); ++ti) {
            if (!bad_set.count(static_cast<int>(ti))) {
                remaining.push_back(mesh.triangles[ti]);
            }
        }
        mesh.triangles.swap(remaining);

        if (polygon.empty()) {
            continue;
        }

        int start = polygon.front().first;
        int current = polygon.front().second;
        std::vector<std::pair<int, int>> ordered;
        ordered.push_back({start, current});
        std::unordered_set<long long> used;
        used.insert((static_cast<long long>(start) << 32) | static_cast<unsigned int>(current));

        for (std::size_t step = 1; step < polygon.size(); ++step) {
            bool found = false;
            for (const auto& edge : polygon) {
                const long long key = (static_cast<long long>(edge.first) << 32) | static_cast<unsigned int>(edge.second);
                if (used.count(key)) {
                    continue;
                }
                if (edge.first == current) {
                    ordered.push_back(edge);
                    current = edge.second;
                    used.insert(key);
                    found = true;
                    break;
                }
                if (edge.second == current) {
                    ordered.push_back({edge.second, edge.first});
                    current = edge.first;
                    used.insert(key);
                    found = true;
                    break;
                }
            }
            if (!found) {
                break;
            }
        }

        for (const auto& edge : ordered) {
            mesh.triangles.push_back({vi, edge.first, edge.second});
        }
    }

    std::vector<std::array<int, 3>> filtered;
    for (const auto& tri : mesh.triangles) {
        if (tri[0] < original_vertex_count || tri[1] < original_vertex_count || tri[2] < original_vertex_count) {
            filtered.push_back(tri);
        }
    }
    mesh.triangles.swap(filtered);
    mesh.vertices.resize(original_vertex_count);
}

int ConstrainedDelaunayTriangulator::findEdgeTriangle(const Triangulation2D& mesh, int v0, int v1) {
    for (std::size_t ti = 0; ti < mesh.triangles.size(); ++ti) {
        const auto& tri = mesh.triangles[ti];
        bool has0 = false;
        bool has1 = false;
        for (int index : tri) {
            if (index == v0) {
                has0 = true;
            }
            if (index == v1) {
                has1 = true;
            }
        }
        if (has0 && has1) {
            return static_cast<int>(ti);
        }
    }
    return -1;
}

int ConstrainedDelaunayTriangulator::oppositeVertex(const std::array<int, 3>& tri, int v0, int v1) {
    for (int index : tri) {
        if (index != v0 && index != v1) {
            return index;
        }
    }
    return -1;
}

bool ConstrainedDelaunayTriangulator::flipEdge(
    Triangulation2D& mesh, int tri0, int tri1, int v0, int v1) {
    const auto& t0 = mesh.triangles[tri0];
    const auto& t1 = mesh.triangles[tri1];
    const int p = oppositeVertex(t0, v0, v1);
    const int q = oppositeVertex(t1, v0, v1);
    if (p < 0 || q < 0) {
        return false;
    }

    const Point2D& P = mesh.vertices[p];
    const Point2D& Q = mesh.vertices[q];
    const Point2D& V0 = mesh.vertices[v0];
    const Point2D& V1 = mesh.vertices[v1];
    if (!inCircumcircle(P, Q, V0, V1) && !inCircumcircle(Q, P, V0, V1)) {
        return false;
    }

    mesh.triangles[tri0] = {p, q, v0};
    mesh.triangles[tri1] = {p, q, v1};
    return true;
}

void ConstrainedDelaunayTriangulator::enforceConstraints(Triangulation2D& mesh) {
    for (const auto& edge : mesh.constrained_edges) {
        int v0 = edge.first;
        int v1 = edge.second;
        for (int iter = 0; iter < 1000; ++iter) {
            int tri_id = findEdgeTriangle(mesh, v0, v1);
            if (tri_id < 0) {
                break;
            }

            const auto& tri = mesh.triangles[tri_id];
            const int opposite = oppositeVertex(tri, v0, v1);
            int neighbor = -1;
            for (std::size_t ti = 0; ti < mesh.triangles.size(); ++ti) {
                if (static_cast<int>(ti) == tri_id) {
                    continue;
                }
                const auto& other = mesh.triangles[ti];
                int shared = 0;
                for (int index : other) {
                    if (index == v0 || index == v1) {
                        ++shared;
                    }
                }
                if (shared == 2 && (other[0] == opposite || other[1] == opposite || other[2] == opposite)) {
                    neighbor = static_cast<int>(ti);
                    break;
                }
            }

            if (neighbor < 0) {
                break;
            }

            if (!flipEdge(mesh, tri_id, neighbor, v0, v1)) {
                break;
            }
        }
    }
}

void ConstrainedDelaunayTriangulator::removeExteriorTriangles(const PLC2D& plc, Triangulation2D& mesh) {
    if (plc.getPolygons().empty()) {
        return;
    }

    const auto& boundary = plc.getPolygons().front().getPoints();
    if (boundary.size() < 3) {
        return;
    }

    std::vector<char> keep(mesh.triangles.size(), 0);
    for (std::size_t ti = 0; ti < mesh.triangles.size(); ++ti) {
        const auto& tri = mesh.triangles[ti];
        const auto& a = mesh.vertices[tri[0]];
        const auto& b = mesh.vertices[tri[1]];
        const auto& c = mesh.vertices[tri[2]];
        const Point2D centroid(
            (a.x_ + b.x_ + c.x_) / 3.0,
            (a.y_ + b.y_ + c.y_) / 3.0);
        if (!pointInPolygon(centroid, boundary)) {
            continue;
        }

        bool contains_hole_seed = false;
        for (const auto& hole : plc.getHoles()) {
            if (pointInTriangle(hole, a, b, c)) {
                contains_hole_seed = true;
                break;
            }
        }
        if (!contains_hole_seed) {
            keep[ti] = 1;
        }
    }

    std::vector<std::array<int, 3>> filtered;
    for (std::size_t ti = 0; ti < mesh.triangles.size(); ++ti) {
        if (keep[ti]) {
            filtered.push_back(mesh.triangles[ti]);
        }
    }
    mesh.triangles.swap(filtered);
}

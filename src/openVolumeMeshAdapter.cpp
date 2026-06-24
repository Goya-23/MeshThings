#include "openVolumeMeshAdapter.h"

#include <OpenVolumeMesh/Geometry/VectorT.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <stdexcept>
#include <vector>

namespace {

using VertexHandle = OpenVolumeMesh::VertexHandle;
using FaceHandle = OpenVolumeMesh::FaceHandle;
using HalfFaceHandle = OpenVolumeMesh::HalfFaceHandle;

std::vector<int> sortedCopy(std::vector<int> vertices) {
    std::sort(vertices.begin(), vertices.end());
    return vertices;
}

class FaceCache {
public:
    explicit FaceCache(VolumeMesh3D& mesh)
        : mesh_(mesh)
    {}

    FaceHandle triangle(int v0, int v1, int v2) {
        const auto key = sortedCopy({v0, v1, v2});
        const auto found = triangles_.find(key);
        if (found != triangles_.end()) {
            return found->second;
        }

        std::vector<VertexHandle> face_vertices = {
            VertexHandle(v0), VertexHandle(v1), VertexHandle(v2)};
        const FaceHandle face = mesh_.add_face(face_vertices);
        triangles_.emplace(std::move(key), face);
        return face;
    }

    FaceHandle quadrilateral(int v0, int v1, int v2, int v3) {
        const auto key = sortedCopy({v0, v1, v2, v3});
        const auto found = quads_.find(key);
        if (found != quads_.end()) {
            return found->second;
        }

        std::vector<VertexHandle> face_vertices = {
            VertexHandle(v0), VertexHandle(v1), VertexHandle(v2), VertexHandle(v3)};
        const FaceHandle face = mesh_.add_face(face_vertices);
        quads_.emplace(std::move(key), face);
        return face;
    }

private:
    VolumeMesh3D& mesh_;
    std::map<std::vector<int>, FaceHandle> triangles_;
    std::map<std::vector<int>, FaceHandle> quads_;
};

HalfFaceHandle orientedHalfFace(VolumeMesh3D& mesh, FaceHandle face, int orientation) {
    return mesh.halfface_handle(face, orientation);
}

double tetVolume6(const std::vector<Point3D>& vertices, const std::array<int, 4>& tet) {
    const Point3D& p0 = vertices[tet[0]];
    const Point3D& p1 = vertices[tet[1]];
    const Point3D& p2 = vertices[tet[2]];
    const Point3D& p3 = vertices[tet[3]];
    return (p1.x_ - p0.x_) * ((p2.y_ - p0.y_) * (p3.z_ - p0.z_) - (p2.z_ - p0.z_) * (p3.y_ - p0.y_)) -
           (p1.y_ - p0.y_) * ((p2.x_ - p0.x_) * (p3.z_ - p0.z_) - (p2.z_ - p0.z_) * (p3.x_ - p0.x_)) +
           (p1.z_ - p0.z_) * ((p2.x_ - p0.x_) * (p3.y_ - p0.y_) - (p2.y_ - p0.y_) * (p3.x_ - p0.x_));
}

void sortRingCounterClockwise(const std::vector<Point3D>& vertices, std::vector<int>& ring) {
    double centroid_x = 0.0;
    double centroid_y = 0.0;
    for (int vertex_id : ring) {
        centroid_x += vertices[vertex_id].x_;
        centroid_y += vertices[vertex_id].y_;
    }
    centroid_x /= static_cast<double>(ring.size());
    centroid_y /= static_cast<double>(ring.size());

    std::sort(ring.begin(), ring.end(), [&](int left_id, int right_id) {
        const Point3D& left = vertices[left_id];
        const Point3D& right = vertices[right_id];
        const double left_angle = std::atan2(left.y_ - centroid_y, left.x_ - centroid_x);
        const double right_angle = std::atan2(right.y_ - centroid_y, right.x_ - centroid_x);
        return left_angle < right_angle;
    });
}

void addWedgeCell(VolumeMesh3D& mesh, FaceCache& faces, const std::array<int, 6>& wedge) {
    const FaceHandle bottom = faces.triangle(wedge[0], wedge[1], wedge[2]);
    const FaceHandle top = faces.triangle(wedge[3], wedge[5], wedge[4]);
    const FaceHandle side0 = faces.quadrilateral(wedge[0], wedge[1], wedge[4], wedge[3]);
    const FaceHandle side1 = faces.quadrilateral(wedge[1], wedge[2], wedge[5], wedge[4]);
    const FaceHandle side2 = faces.quadrilateral(wedge[2], wedge[0], wedge[3], wedge[5]);

    std::vector<HalfFaceHandle> halffaces = {
        orientedHalfFace(mesh, bottom, 0),
        orientedHalfFace(mesh, top, 1),
        orientedHalfFace(mesh, side0, 1),
        orientedHalfFace(mesh, side1, 1),
        orientedHalfFace(mesh, side2, 1),
    };
    mesh.add_cell(halffaces, false);
}

}  // namespace

SurfaceMesh2D OpenVolumeMeshAdapter::buildSurfaceMesh(const Triangulation2D& triangulation) {
    SurfaceMesh2D mesh;

    std::vector<VertexHandle> handles;
    handles.reserve(triangulation.vertices.size());
    for (const auto& vertex : triangulation.vertices) {
        handles.push_back(mesh.add_vertex(OpenVolumeMesh::Geometry::Vec2d(vertex.x_, vertex.y_)));
    }

    for (const auto& triangle : triangulation.triangles) {
        std::vector<VertexHandle> face_vertices = {
            handles[triangle[0]], handles[triangle[1]], handles[triangle[2]]};
        mesh.add_face(face_vertices);
    }

    return mesh;
}

Triangulation2D OpenVolumeMeshAdapter::extractTriangulation(const SurfaceMesh2D& mesh) {
    Triangulation2D triangulation;
    triangulation.vertices.reserve(mesh.n_vertices());
    for (OpenVolumeMesh::VertexIter vertex_it = mesh.vertices_begin();
         vertex_it != mesh.vertices_end();
         ++vertex_it) {
        const auto position = mesh.vertex(*vertex_it);
        triangulation.vertices.emplace_back(position[0], position[1]);
    }

    for (OpenVolumeMesh::FaceIter face_it = mesh.faces_begin(); face_it != mesh.faces_end(); ++face_it) {
        std::array<int, 3> triangle{};
        int corner = 0;
        for (OpenVolumeMesh::FaceVertexIter fv_it = mesh.fv_iter(*face_it); fv_it.valid(); ++fv_it) {
            triangle[corner++] = fv_it->idx();
        }
        if (corner == 3) {
            triangulation.triangles.push_back(triangle);
        }
    }

    return triangulation;
}

VolumeMesh3D OpenVolumeMeshAdapter::buildWedgeVolumeMesh(const WedgeMesh3D& wedge_mesh) {
    VolumeMesh3D mesh;
    for (const auto& vertex : wedge_mesh.vertices) {
        mesh.add_vertex(OpenVolumeMesh::Geometry::Vec3d(vertex.x_, vertex.y_, vertex.z_));
    }

    FaceCache faces(mesh);
    for (const auto& wedge : wedge_mesh.wedges) {
        addWedgeCell(mesh, faces, wedge);
    }
    return mesh;
}

WedgeMesh3D OpenVolumeMeshAdapter::extractWedgeMesh(const VolumeMesh3D& mesh, double conductivity) {
    WedgeMesh3D wedge_mesh;
    wedge_mesh.conductivity = conductivity;
    wedge_mesh.vertices.reserve(mesh.n_vertices());
    for (OpenVolumeMesh::VertexIter vertex_it = mesh.vertices_begin();
         vertex_it != mesh.vertices_end();
         ++vertex_it) {
        const auto position = mesh.vertex(*vertex_it);
        wedge_mesh.vertices.emplace_back(position[0], position[1], position[2]);
    }

    for (OpenVolumeMesh::CellIter cell_it = mesh.cells_begin(); cell_it != mesh.cells_end(); ++cell_it) {
        std::vector<int> cell_vertices;
        for (OpenVolumeMesh::CellVertexIter cv_it = mesh.cv_iter(*cell_it); cv_it.valid(); ++cv_it) {
            cell_vertices.push_back(cv_it->idx());
        }
        if (cell_vertices.size() != 6) {
            continue;
        }

        std::array<int, 6> wedge{};
        double min_z = wedge_mesh.vertices[cell_vertices.front()].z_;
        double max_z = min_z;
        for (int vertex_id : cell_vertices) {
            min_z = std::min(min_z, wedge_mesh.vertices[vertex_id].z_);
            max_z = std::max(max_z, wedge_mesh.vertices[vertex_id].z_);
        }

        std::vector<int> bottom;
        std::vector<int> top;
        for (int vertex_id : cell_vertices) {
            if (std::abs(wedge_mesh.vertices[vertex_id].z_ - min_z) < 1.0e-10) {
                bottom.push_back(vertex_id);
            } else {
                top.push_back(vertex_id);
            }
        }
        if (bottom.size() != 3 || top.size() != 3) {
            continue;
        }

        sortRingCounterClockwise(wedge_mesh.vertices, bottom);

        const auto matchTopVertex = [&](int bottom_id) {
            const Point3D& bottom_point = wedge_mesh.vertices[bottom_id];
            for (int top_id : top) {
                const Point3D& top_point = wedge_mesh.vertices[top_id];
                if (std::abs(bottom_point.x_ - top_point.x_) < 1.0e-10 &&
                    std::abs(bottom_point.y_ - top_point.y_) < 1.0e-10) {
                    return top_id;
                }
            }
            throw std::runtime_error("Failed to match wedge top vertex to bottom footprint");
        };

        wedge[0] = bottom[0];
        wedge[1] = bottom[1];
        wedge[2] = bottom[2];
        wedge[3] = matchTopVertex(bottom[0]);
        wedge[4] = matchTopVertex(bottom[1]);
        wedge[5] = matchTopVertex(bottom[2]);

        const std::array<int, 4> probe_tet = {wedge[0], wedge[1], wedge[2], wedge[5]};
        if (tetVolume6(wedge_mesh.vertices, probe_tet) < 0.0) {
            std::swap(wedge[1], wedge[2]);
            std::swap(wedge[4], wedge[5]);
        }
        wedge_mesh.wedges.push_back(wedge);
    }

    if (wedge_mesh.wedgeCount() == 0) {
        throw std::runtime_error("OpenVolumeMesh volume does not contain wedge cells");
    }
    return wedge_mesh;
}

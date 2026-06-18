#include "openVolumeMeshAdapter.h"

#include <OpenVolumeMesh/Geometry/VectorT.hh>

SurfaceMesh2D OpenVolumeMeshAdapter::buildSurfaceMesh(const Triangulation2D& triangulation) {
    SurfaceMesh2D mesh;

    std::vector<OpenVolumeMesh::VertexHandle> handles;
    handles.reserve(triangulation.vertices.size());
    for (const auto& vertex : triangulation.vertices) {
        handles.push_back(mesh.add_vertex(OpenVolumeMesh::Geometry::Vec2d(vertex.x_, vertex.y_)));
    }

    for (const auto& triangle : triangulation.triangles) {
        std::vector<OpenVolumeMesh::VertexHandle> face_vertices = {
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

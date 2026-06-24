#include "wedgeMesh.h"

#include <stdexcept>

WedgeMesh3D WedgeExtruder::extrude(
    const Triangulation2D& base, int nz, double height, double conductivity) {
    if (nz < 1) {
        throw std::runtime_error("Wedge extrusion requires at least one layer (nz >= 1)");
    }
    if (height <= 0.0) {
        throw std::runtime_error("Wedge extrusion height must be positive");
    }

    WedgeMesh3D mesh;
    mesh.conductivity = conductivity;

    const int vertices_per_layer = static_cast<int>(base.vertexCount());
    const int layers = nz + 1;
    mesh.vertices.reserve(static_cast<std::size_t>(vertices_per_layer * layers));

    for (int layer = 0; layer <= nz; ++layer) {
        const double z = height * static_cast<double>(layer) / static_cast<double>(nz);
        for (const auto& vertex : base.vertices) {
            mesh.vertices.emplace_back(vertex.x_, vertex.y_, z);
        }
    }

    auto layer_vertex = [vertices_per_layer](int layer, int base_vertex) {
        return layer * vertices_per_layer + base_vertex;
    };

    mesh.wedges.reserve(base.triangleCount() * static_cast<std::size_t>(nz));
    for (int layer = 0; layer < nz; ++layer) {
        for (const auto& triangle : base.triangles) {
            mesh.wedges.push_back({
                layer_vertex(layer, triangle[0]),
                layer_vertex(layer, triangle[1]),
                layer_vertex(layer, triangle[2]),
                layer_vertex(layer + 1, triangle[0]),
                layer_vertex(layer + 1, triangle[1]),
                layer_vertex(layer + 1, triangle[2]),
            });
        }
    }

    return mesh;
}

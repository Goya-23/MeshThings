#include "wedgeMesh.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace {

constexpr double kGeometryTolerance = 1.0e-10;

bool pointOnSegment(const Point2D& point, const Point2D& a, const Point2D& b) {
    const double cross = (point.x_ - a.x_) * (b.y_ - a.y_) - (point.y_ - a.y_) * (b.x_ - a.x_);
    if (std::abs(cross) > kGeometryTolerance) {
        return false;
    }

    const double xmin = std::min(a.x_, b.x_) - kGeometryTolerance;
    const double xmax = std::max(a.x_, b.x_) + kGeometryTolerance;
    const double ymin = std::min(a.y_, b.y_) - kGeometryTolerance;
    const double ymax = std::max(a.y_, b.y_) + kGeometryTolerance;
    return point.x_ >= xmin && point.x_ <= xmax && point.y_ >= ymin && point.y_ <= ymax;
}

std::vector<char> findBaseBoundaryVertices(const Triangulation2D& base) {
    std::vector<char> boundary(base.vertexCount(), 0);
    for (const auto& edge : base.constrained_edges) {
        if (edge.first < 0 || edge.first >= static_cast<int>(boundary.size()) ||
            edge.second < 0 || edge.second >= static_cast<int>(boundary.size())) {
            continue;
        }

        boundary[edge.first] = 1;
        boundary[edge.second] = 1;
        const Point2D& a = base.vertices[edge.first];
        const Point2D& b = base.vertices[edge.second];
        for (std::size_t vertex_id = 0; vertex_id < base.vertices.size(); ++vertex_id) {
            if (pointOnSegment(base.vertices[vertex_id], a, b)) {
                boundary[vertex_id] = 1;
            }
        }
    }
    return boundary;
}

bool contains(const MaterialRegion2D& region, double x, double y, double z) {
    return x >= region.xmin - kGeometryTolerance && x <= region.xmax + kGeometryTolerance &&
           y >= region.ymin - kGeometryTolerance && y <= region.ymax + kGeometryTolerance &&
           z >= region.zmin - kGeometryTolerance && z <= region.zmax + kGeometryTolerance;
}

}  // namespace

WedgeMesh3D WedgeExtruder::extrude(
    const Triangulation2D& base,
    int nz,
    double height,
    double conductivity,
    const std::vector<MaterialRegion2D>& material_regions) {
    if (nz < 1) {
        throw std::runtime_error("Wedge extrusion requires at least one layer (nz >= 1)");
    }
    if (height <= 0.0) {
        throw std::runtime_error("Wedge extrusion height must be positive");
    }

    WedgeMesh3D mesh;
    mesh.conductivity = conductivity;
    mesh.material_names.push_back("default");
    for (const auto& material : material_regions) {
        mesh.material_names.push_back(material.name);
    }

    const int vertices_per_layer = static_cast<int>(base.vertexCount());
    const int layers = nz + 1;
    mesh.vertices.reserve(static_cast<std::size_t>(vertices_per_layer * layers));
    mesh.boundary_vertices.assign(static_cast<std::size_t>(vertices_per_layer * layers), 0);
    const std::vector<char> base_boundary = findBaseBoundaryVertices(base);

    for (int layer = 0; layer <= nz; ++layer) {
        const double z = height * static_cast<double>(layer) / static_cast<double>(nz);
        for (std::size_t base_vertex = 0; base_vertex < base.vertices.size(); ++base_vertex) {
            const auto& vertex = base.vertices[base_vertex];
            mesh.vertices.emplace_back(vertex.x_, vertex.y_, z);
            const int vertex_id = layer * vertices_per_layer + static_cast<int>(base_vertex);
            if (layer == 0 || layer == nz || (base_vertex < base_boundary.size() && base_boundary[base_vertex])) {
                mesh.boundary_vertices[vertex_id] = 1;
            }
        }
    }

    auto layer_vertex = [vertices_per_layer](int layer, int base_vertex) {
        return layer * vertices_per_layer + base_vertex;
    };

    mesh.wedges.reserve(base.triangleCount() * static_cast<std::size_t>(nz));
    mesh.wedge_conductivity.reserve(base.triangleCount() * static_cast<std::size_t>(nz));
    mesh.wedge_material_id.reserve(base.triangleCount() * static_cast<std::size_t>(nz));
    for (int layer = 0; layer < nz; ++layer) {
        for (const auto& triangle : base.triangles) {
            const auto wedge = std::array<int, 6>{
                layer_vertex(layer, triangle[0]),
                layer_vertex(layer, triangle[1]),
                layer_vertex(layer, triangle[2]),
                layer_vertex(layer + 1, triangle[0]),
                layer_vertex(layer + 1, triangle[1]),
                layer_vertex(layer + 1, triangle[2]),
            };
            mesh.wedges.push_back({
                wedge[0],
                wedge[1],
                wedge[2],
                wedge[3],
                wedge[4],
                wedge[5],
            });

            double centroid_x = 0.0;
            double centroid_y = 0.0;
            double centroid_z = 0.0;
            for (int vertex_id : wedge) {
                centroid_x += mesh.vertices[vertex_id].x_;
                centroid_y += mesh.vertices[vertex_id].y_;
                centroid_z += mesh.vertices[vertex_id].z_;
            }
            centroid_x /= 6.0;
            centroid_y /= 6.0;
            centroid_z /= 6.0;

            double wedge_conductivity = conductivity;
            int material_id = 0;
            for (std::size_t region_id = 0; region_id < material_regions.size(); ++region_id) {
                if (contains(material_regions[region_id], centroid_x, centroid_y, centroid_z)) {
                    wedge_conductivity = material_regions[region_id].conductivity;
                    material_id = static_cast<int>(region_id) + 1;
                }
            }
            mesh.wedge_conductivity.push_back(wedge_conductivity);
            mesh.wedge_material_id.push_back(material_id);
        }
    }

    return mesh;
}

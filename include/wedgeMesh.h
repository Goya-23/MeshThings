#ifndef _WEDGE_MESH_H_
#define _WEDGE_MESH_H_

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "point.h"
#include "plc.h"
#include "triangulation2d.h"

struct WedgeMesh3D {
    std::vector<Point3D> vertices;
    std::vector<std::array<int, 6>> wedges;
    std::vector<double> wedge_conductivity;
    std::vector<int> wedge_material_id;
    std::vector<std::string> material_names;
    std::vector<char> boundary_vertices;
    double conductivity = 1.0;

    std::size_t vertexCount() const { return vertices.size(); }
    std::size_t wedgeCount() const { return wedges.size(); }
    double conductivityForWedge(std::size_t wedge_id) const {
        if (wedge_id < wedge_conductivity.size()) {
            return wedge_conductivity[wedge_id];
        }
        return conductivity;
    }
};

class WedgeExtruder {
public:
    static WedgeMesh3D extrude(
        const Triangulation2D& base,
        int nz,
        double height,
        double conductivity = 1.0,
        const std::vector<MaterialRegion2D>& material_regions = {});
};

#endif

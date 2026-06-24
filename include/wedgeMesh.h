#ifndef _WEDGE_MESH_H_
#define _WEDGE_MESH_H_

#include <array>
#include <cstddef>
#include <vector>

#include "point.h"
#include "triangulation2d.h"

struct WedgeMesh3D {
    std::vector<Point3D> vertices;
    std::vector<std::array<int, 6>> wedges;
    double conductivity = 1.0;

    std::size_t vertexCount() const { return vertices.size(); }
    std::size_t wedgeCount() const { return wedges.size(); }
};

class WedgeExtruder {
public:
    static WedgeMesh3D extrude(const Triangulation2D& base, int nz, double height, double conductivity = 1.0);
};

#endif

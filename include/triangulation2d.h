#ifndef _TRIANGULATION_2D_H_
#define _TRIANGULATION_2D_H_

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#include "point.h"

struct Triangulation2D {
    std::vector<Point2D> vertices;
    std::vector<std::array<int, 3>> triangles;
    std::vector<std::pair<int, int>> constrained_edges;

    std::size_t vertexCount() const { return vertices.size(); }
    std::size_t triangleCount() const { return triangles.size(); }
};

#endif

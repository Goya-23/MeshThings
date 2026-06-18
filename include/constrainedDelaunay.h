#ifndef _CONSTRAINED_DELAUNAY_H_
#define _CONSTRAINED_DELAUNAY_H_

#include <memory>
#include <vector>

#include "plc.h"
#include "triangulation2d.h"

class ConstrainedDelaunayTriangulator {
public:
    Triangulation2D triangulate(const PLC2D& plc, int nx, int ny);

private:
    void addInteriorPoints(const PLC2D& plc, int nx, int ny, Triangulation2D& mesh);
    void buildBoundaryConstraints(const PLC2D& plc, Triangulation2D& mesh);
    void bowyerWatson(Triangulation2D& mesh);
    void enforceConstraints(Triangulation2D& mesh);
    void removeExteriorTriangles(const PLC2D& plc, Triangulation2D& mesh);

    static bool inCircumcircle(const Point2D& p, const Point2D& a, const Point2D& b, const Point2D& c);
    static int findEdgeTriangle(const Triangulation2D& mesh, int v0, int v1);
    static int oppositeVertex(const std::array<int, 3>& tri, int v0, int v1);
    static bool flipEdge(Triangulation2D& mesh, int tri0, int tri1, int v0, int v1);
};

#endif

#include "concreteDelaunayTriangulation.h"

#include "constrainedDelaunay.h"

void IncrementalDelaunayTriangulation::triangulate(std::shared_ptr<PLC2D> plc) {
    if (!plc) {
        return;
    }
    ConstrainedDelaunayTriangulator triangulator;
    result_ = triangulator.triangulate(*plc, nx_, ny_);
}

void SweepLineDelaunayTriangulation::triangulate(std::shared_ptr<PLC2D> plc) {
    IncrementalDelaunayTriangulation fallback;
    fallback.setGridDensity(nx_, ny_);
    fallback.triangulate(plc);
    result_ = fallback.result();
}

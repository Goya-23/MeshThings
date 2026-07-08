#ifndef _CONCRETE_DELAUNAY_TRIANGULATION_H_
#define _CONCRETE_DELAUNAY_TRIANGULATION_H_
#include <memory>
#include "plc.h"
#include "delaunayTriangulation.h"

class IncrementalDelaunayTriangulation : public DelaunayTriangulation {
public:
    IncrementalDelaunayTriangulation() {}
    virtual ~IncrementalDelaunayTriangulation() {}

    void setGridDensity(int nx, int ny) {
        nx_ = nx;
        ny_ = ny;
    }

    void triangulate(std::shared_ptr<PLC2D> plc) override;

private:
    int nx_ = 8;
    int ny_ = 8;
};

class SweepLineDelaunayTriangulation : public DelaunayTriangulation {
public:
    SweepLineDelaunayTriangulation() {}
    virtual ~SweepLineDelaunayTriangulation() {}

    void setGridDensity(int nx, int ny) {
        nx_ = nx;
        ny_ = ny;
    }

    void triangulate(std::shared_ptr<PLC2D> plc) override;

private:
    int nx_ = 8;
    int ny_ = 8;
};

#endif  // _CONCRETE_DELAUNAY_TRIANGULATION_H_

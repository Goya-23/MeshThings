#ifndef _DELAUNAY_TRIANGULATION_H_
#define _DELAUNAY_TRIANGULATION_H_
#include <memory>
#include "plc.h"
#include "triangulation2d.h"

class DelaunayTriangulation {
public:
    DelaunayTriangulation() {}
    virtual ~DelaunayTriangulation() {}

    virtual void triangulate(std::shared_ptr<PLC2D> plc) = 0;
    const Triangulation2D& result() const { return result_; }

protected:
    Triangulation2D result_;
};

#endif  // _DELAUNAY_TRIANGULATION_H_

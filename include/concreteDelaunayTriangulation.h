#ifndef _CONCRETE_DELAUNAY_TRIANGULATION_H_
#define _CONCRETE_DELAUNAY_TRIANGULATION_H_
#include <memory>
#include "plc.h"



class IncrementalDelaunayTriangulation
{
public:
    IncrementalDelaunayTriangulation() {}
    virtual ~IncrementalDelaunayTriangulation() {}
    
    void triangulate(std::shared_ptr<PLC2D> plc);



};



class SweepLineDelaunayTriangulation
{
public:
    SweepLineDelaunayTriangulation() {}
    virtual ~SweepLineDelaunayTriangulation() {}

    void triangulate(std::shared_ptr<PLC2D> plc);
};

#endif // _CONCRETE_DELAUNAY_TRIANGULATION_H_
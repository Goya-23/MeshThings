#ifndef _DELAUNAY_TRIANGULATION_H_
#define _DELAUNAY_TRIANGULATION_H_
#include <memory>
#include "plc.h"



class DelaunayTriangulation
{
public:
    DelaunayTriangulation() {}
    virtual ~DelaunayTriangulation() {}
    
    virtual void triangulate(std::shared_ptr<PLC2D> plc) = 0;



};





#endif // _DELAUNAY_TRIANGULATION_H_
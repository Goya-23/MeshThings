#ifndef _PLC_H_
#define _PLC_H_
#include <vector>
#include <memory>
#include "polygon.h"

class PLC2D
{
public:
    PLC2D(const std::vector<Polygon>& polygons, const std::vector<Point2D>& holes) 
        : m_polygons(polygons), m_holes(holes)
    {}

    const std::vector<Polygon>& getPolygons() const { return m_polygons; }
    const std::vector<Point2D>& getHoles() const { return m_holes; }

private:
    std::vector<Polygon> m_polygons;
    std::vector<Point2D> m_holes;

}



#endif // _PLC_H_
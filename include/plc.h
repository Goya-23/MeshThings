#ifndef _PLC_H_
#define _PLC_H_
#include <vector>
#include <memory>
#include "polygon.h"

class PLC2D {
public:
    PLC2D(const std::vector<Polygon>& polygons, const std::vector<Point2D>& holes, double extrusion_height = 1.0)
        : m_polygons(polygons), m_holes(holes), m_extrusion_height(extrusion_height)
    {}

    const std::vector<Polygon>& getPolygons() const { return m_polygons; }
    const std::vector<Point2D>& getHoles() const { return m_holes; }
    double getExtrusionHeight() const { return m_extrusion_height; }

private:
    std::vector<Polygon> m_polygons;
    std::vector<Point2D> m_holes;
    double m_extrusion_height;
};

#endif  // _PLC_H_

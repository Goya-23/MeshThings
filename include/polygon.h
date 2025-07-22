#ifndef _POLYGON_H_
#define _POLYGON_H_

#include "point.h"

class Polygon 
{   
public:
    Polygon(const std::vector<Point2D>& points)
        : m_points(points)    
    {}
    ~Polygon() {}

    const std::vector<Point2D>& getPoints() const { return m_points; }

private:
    std::vector<Point2D> m_points;
};





#endif // _POLYGON_H_
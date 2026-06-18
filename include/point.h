#ifndef _POINT_H_
#define _POINT_H_

class Point2D 
{
public:
    Point2D()
        : x_(0.0), y_(0.0)
    {}

    Point2D(double x, double y)
        : x_(x), y_(y)
    {}

    double x_;
    double y_;
};



#endif // _POINT_H_
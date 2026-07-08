#ifndef _POINT_H_
#define _POINT_H_

class Point2D {
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

class Point3D {
public:
    Point3D()
        : x_(0.0), y_(0.0), z_(0.0)
    {}

    Point3D(double x, double y, double z)
        : x_(x), y_(y), z_(z)
    {}

    Point3D(const Point2D& point, double z)
        : x_(point.x_), y_(point.y_), z_(z)
    {}

    double x_;
    double y_;
    double z_;
};

#endif  // _POINT_H_

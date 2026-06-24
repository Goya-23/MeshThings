#include "plcParser.h"

#include <fstream>
#include <sstream>
#include <stdexcept>

std::shared_ptr<PLC2D> PLCParser::parse(const char* path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error(std::string("Unable to open PLC file: ") + path);
    }

    std::vector<Point2D> boundary;
    int boundary_count = 0;
    input >> boundary_count;
    boundary.reserve(static_cast<std::size_t>(boundary_count));
    for (int i = 0; i < boundary_count; ++i) {
        double x = 0.0;
        double y = 0.0;
        input >> x >> y;
        boundary.emplace_back(x, y);
    }

    int hole_count = 0;
    input >> hole_count;
    std::vector<Point2D> holes;
    holes.reserve(static_cast<std::size_t>(hole_count));
    for (int i = 0; i < hole_count; ++i) {
        double x = 0.0;
        double y = 0.0;
        input >> x >> y;
        holes.emplace_back(x, y);
    }

    double extrusion_height = 1.0;
    input >> extrusion_height;

    if (boundary.size() < 3) {
        throw std::runtime_error("PLC boundary must contain at least three vertices");
    }

    std::vector<Polygon> polygons;
    polygons.emplace_back(boundary);
    return std::make_shared<PLC2D>(polygons, holes, extrusion_height);
}

std::shared_ptr<PLC2D> PLCParser::makeUnitSquare() {
    std::vector<Point2D> boundary = {
        Point2D(0.0, 0.0),
        Point2D(1.0, 0.0),
        Point2D(1.0, 1.0),
        Point2D(0.0, 1.0),
    };
    std::vector<Polygon> polygons;
    polygons.emplace_back(boundary);
    return std::make_shared<PLC2D>(polygons, std::vector<Point2D>{}, 1.0);
}

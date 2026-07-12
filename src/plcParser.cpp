#include "plcParser.h"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

class TokenReader {
public:
    explicit TokenReader(std::istream& input)
        : input_(input)
    {}

    bool next(std::string& token) {
        while (input_ >> token) {
            if (!token.empty() && token[0] == '#') {
                std::string ignored;
                std::getline(input_, ignored);
                continue;
            }
            return true;
        }
        return false;
    }

    template <typename T>
    bool nextValue(T& value) {
        std::string token;
        if (!next(token)) {
            return false;
        }
        std::istringstream stream(token);
        stream >> value;
        return !stream.fail();
    }

private:
    std::istream& input_;
};

}  // namespace

std::shared_ptr<PLC2D> PLCParser::parse(const char* path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error(std::string("Unable to open PLC file: ") + path);
    }
    TokenReader reader(input);

    std::vector<Point2D> boundary;
    int boundary_count = 0;
    if (!reader.nextValue(boundary_count)) {
        throw std::runtime_error("PLC file is missing boundary vertex count");
    }
    if (boundary_count < 3) {
        throw std::runtime_error("PLC boundary must contain at least three vertices");
    }
    boundary.reserve(static_cast<std::size_t>(boundary_count));
    for (int i = 0; i < boundary_count; ++i) {
        double x = 0.0;
        double y = 0.0;
        if (!reader.nextValue(x) || !reader.nextValue(y)) {
            throw std::runtime_error("PLC file ended while reading boundary vertices");
        }
        boundary.emplace_back(x, y);
    }

    int hole_count = 0;
    if (!reader.nextValue(hole_count)) {
        throw std::runtime_error("PLC file is missing hole seed count");
    }
    if (hole_count < 0) {
        throw std::runtime_error("PLC hole seed count cannot be negative");
    }
    std::vector<Point2D> holes;
    holes.reserve(static_cast<std::size_t>(hole_count));
    for (int i = 0; i < hole_count; ++i) {
        double x = 0.0;
        double y = 0.0;
        if (!reader.nextValue(x) || !reader.nextValue(y)) {
            throw std::runtime_error("PLC file ended while reading hole seed points");
        }
        holes.emplace_back(x, y);
    }

    double extrusion_height = 1.0;
    reader.nextValue(extrusion_height);

    double default_conductivity = 1.0;
    if (!reader.nextValue(default_conductivity)) {
        default_conductivity = 1.0;
    }

    int material_region_count = 0;
    if (!reader.nextValue(material_region_count)) {
        material_region_count = 0;
    }
    if (material_region_count < 0) {
        throw std::runtime_error("PLC material region count cannot be negative");
    }
    std::vector<MaterialRegion2D> material_regions;
    material_regions.reserve(static_cast<std::size_t>(material_region_count));
    for (int i = 0; i < material_region_count; ++i) {
        MaterialRegion2D region;
        if (!reader.next(region.name) ||
            !reader.nextValue(region.conductivity) ||
            !reader.nextValue(region.xmin) ||
            !reader.nextValue(region.xmax) ||
            !reader.nextValue(region.ymin) ||
            !reader.nextValue(region.ymax) ||
            !reader.nextValue(region.zmin) ||
            !reader.nextValue(region.zmax)) {
            throw std::runtime_error("PLC file ended while reading material regions");
        }
        if (region.xmin > region.xmax) {
            std::swap(region.xmin, region.xmax);
        }
        if (region.ymin > region.ymax) {
            std::swap(region.ymin, region.ymax);
        }
        if (region.zmin > region.zmax) {
            std::swap(region.zmin, region.zmax);
        }
        material_regions.push_back(region);
    }

    std::vector<Polygon> polygons;
    polygons.emplace_back(boundary);
    return std::make_shared<PLC2D>(
        polygons,
        holes,
        extrusion_height,
        default_conductivity,
        material_regions);
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
    return std::make_shared<PLC2D>(polygons, std::vector<Point2D>{}, 1.0, 1.0);
}

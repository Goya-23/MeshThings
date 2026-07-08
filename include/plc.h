#ifndef _PLC_H_
#define _PLC_H_
#include <string>
#include <vector>
#include <memory>
#include "polygon.h"

struct MaterialRegion2D {
    std::string name;
    double conductivity = 1.0;
    double xmin = 0.0;
    double xmax = 0.0;
    double ymin = 0.0;
    double ymax = 0.0;
    double zmin = 0.0;
    double zmax = 0.0;
};

class PLC2D {
public:
    PLC2D(const std::vector<Polygon>& polygons,
          const std::vector<Point2D>& holes,
          double extrusion_height = 1.0,
          double default_conductivity = 1.0,
          const std::vector<MaterialRegion2D>& material_regions = {})
        : m_polygons(polygons),
          m_holes(holes),
          m_extrusion_height(extrusion_height),
          m_default_conductivity(default_conductivity),
          m_material_regions(material_regions)
    {}

    const std::vector<Polygon>& getPolygons() const { return m_polygons; }
    const std::vector<Point2D>& getHoles() const { return m_holes; }
    double getExtrusionHeight() const { return m_extrusion_height; }
    double getDefaultConductivity() const { return m_default_conductivity; }
    const std::vector<MaterialRegion2D>& getMaterialRegions() const { return m_material_regions; }

private:
    std::vector<Polygon> m_polygons;
    std::vector<Point2D> m_holes;
    double m_extrusion_height;
    double m_default_conductivity;
    std::vector<MaterialRegion2D> m_material_regions;
};

#endif  // _PLC_H_

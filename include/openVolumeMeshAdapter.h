#ifndef _OPEN_VOLUME_MESH_ADAPTER_H_
#define _OPEN_VOLUME_MESH_ADAPTER_H_

#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>

#include <string>
#include <vector>

#include "triangulation2d.h"
#include "wedgeMesh.h"

using SurfaceMesh2D = OpenVolumeMesh::GeometricPolyhedralMeshV2d;
using VolumeMesh3D = OpenVolumeMesh::GeometricPolyhedralMeshV3d;

class OpenVolumeMeshAdapter {
public:
    static SurfaceMesh2D buildSurfaceMesh(const Triangulation2D& triangulation);
    static Triangulation2D extractTriangulation(const SurfaceMesh2D& mesh);

    static VolumeMesh3D buildWedgeVolumeMesh(const WedgeMesh3D& wedge_mesh);
    static WedgeMesh3D extractWedgeMesh(
        const VolumeMesh3D& mesh,
        double conductivity = 1.0,
        const std::vector<double>& wedge_conductivity = {},
        const std::vector<int>& wedge_material_id = {},
        const std::vector<std::string>& material_names = {},
        const std::vector<char>& boundary_vertices = {});
};

#endif

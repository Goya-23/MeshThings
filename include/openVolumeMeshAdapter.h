#ifndef _OPEN_VOLUME_MESH_ADAPTER_H_
#define _OPEN_VOLUME_MESH_ADAPTER_H_

#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>

#include "triangulation2d.h"

using SurfaceMesh2D = OpenVolumeMesh::GeometricPolyhedralMeshV2d;

class OpenVolumeMeshAdapter {
public:
    static SurfaceMesh2D buildSurfaceMesh(const Triangulation2D& triangulation);
    static Triangulation2D extractTriangulation(const SurfaceMesh2D& mesh);
};

#endif

#ifndef _MESH_PARTITIONER_H_
#define _MESH_PARTITIONER_H_

#include <vector>

#include "triangulation2d.h"

struct MeshPartition {
    int num_parts = 0;
    std::vector<int> element_part;
    std::vector<int> vertex_part;
};

class MeshPartitioner {
public:
    static MeshPartition partitionNodal(const Triangulation2D& mesh, int num_parts);
    static Triangulation2D reorderByPartition(Triangulation2D mesh, MeshPartition& partition);
};

#endif

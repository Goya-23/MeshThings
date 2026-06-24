#include "meshPartitioner.h"

#include <metis.h>

#include <algorithm>
#include <stdexcept>

MeshPartition MeshPartitioner::partitionNodal(const Triangulation2D& mesh, int num_parts) {
    if (num_parts < 1) {
        throw std::runtime_error("METIS partition requires at least one part");
    }

    MeshPartition partition;
    partition.num_parts = num_parts;
    partition.element_part.resize(mesh.triangleCount());
    partition.vertex_part.resize(mesh.vertexCount(), 0);

    if (mesh.triangleCount() == 0) {
        return partition;
    }

    if (num_parts == 1) {
        std::fill(partition.element_part.begin(), partition.element_part.end(), 0);
        std::fill(partition.vertex_part.begin(), partition.vertex_part.end(), 0);
        return partition;
    }

    std::vector<idx_t> eptr(mesh.triangleCount() + 1);
    std::vector<idx_t> eind(mesh.triangleCount() * 3);
    for (std::size_t element = 0; element < mesh.triangleCount(); ++element) {
        eptr[element] = static_cast<idx_t>(3 * element);
        for (int corner = 0; corner < 3; ++corner) {
            eind[3 * element + corner] = static_cast<idx_t>(mesh.triangles[element][corner]);
        }
    }
    eptr[mesh.triangleCount()] = static_cast<idx_t>(eind.size());

    idx_t ne = static_cast<idx_t>(mesh.triangleCount());
    idx_t nn = static_cast<idx_t>(mesh.vertexCount());
    idx_t ncommon = 2;
    idx_t nparts = static_cast<idx_t>(num_parts);
    idx_t objval = 0;

    std::vector<idx_t> epart(mesh.triangleCount());
    std::vector<idx_t> npart(mesh.vertexCount());
    std::vector<idx_t> options(METIS_NOPTIONS);
    METIS_SetDefaultOptions(options.data());
    options[METIS_OPTION_NUMBERING] = 0;

    const int status = METIS_PartMeshNodal(
        &ne,
        &nn,
        eptr.data(),
        eind.data(),
        nullptr,
        nullptr,
        &nparts,
        nullptr,
        options.data(),
        &objval,
        epart.data(),
        npart.data());

    if (status != METIS_OK) {
        throw std::runtime_error("METIS_PartMeshNodal failed");
    }

    for (std::size_t i = 0; i < epart.size(); ++i) {
        partition.element_part[i] = static_cast<int>(epart[i]);
    }
    for (std::size_t i = 0; i < npart.size(); ++i) {
        partition.vertex_part[i] = static_cast<int>(npart[i]);
    }

    return partition;
}

Triangulation2D MeshPartitioner::reorderByPartition(Triangulation2D mesh, MeshPartition& partition) {
    const int vertex_count = static_cast<int>(mesh.vertexCount());
    std::vector<int> order(vertex_count);
    for (int i = 0; i < vertex_count; ++i) {
        order[i] = i;
    }

    std::stable_sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        if (partition.vertex_part[lhs] != partition.vertex_part[rhs]) {
            return partition.vertex_part[lhs] < partition.vertex_part[rhs];
        }
        return lhs < rhs;
    });

    std::vector<int> old_to_new(vertex_count);
    std::vector<int> new_part(vertex_count);
    for (int new_index = 0; new_index < vertex_count; ++new_index) {
        const int old_index = order[new_index];
        old_to_new[old_index] = new_index;
        new_part[new_index] = partition.vertex_part[old_index];
    }

    Triangulation2D reordered;
    reordered.vertices.resize(mesh.vertices.size());
    for (int old_index = 0; old_index < vertex_count; ++old_index) {
        reordered.vertices[old_to_new[old_index]] = mesh.vertices[old_index];
    }

    reordered.triangles.reserve(mesh.triangles.size());
    for (const auto& triangle : mesh.triangles) {
        reordered.triangles.push_back({
            old_to_new[triangle[0]],
            old_to_new[triangle[1]],
            old_to_new[triangle[2]],
        });
    }

    reordered.constrained_edges.reserve(mesh.constrained_edges.size());
    for (const auto& edge : mesh.constrained_edges) {
        reordered.constrained_edges.emplace_back(old_to_new[edge.first], old_to_new[edge.second]);
    }

    partition.vertex_part = std::move(new_part);
    partition.element_part.resize(mesh.triangleCount());
    for (std::size_t element = 0; element < mesh.triangleCount(); ++element) {
        const auto& tri = mesh.triangles[element];
        partition.element_part[element] = partition.vertex_part[old_to_new[tri[0]]];
    }

    return reordered;
}

MeshPartition MeshPartitioner::partitionNodal(const WedgeMesh3D& mesh, int num_parts) {
    if (num_parts < 1) {
        throw std::runtime_error("METIS partition requires at least one part");
    }

    MeshPartition partition;
    partition.num_parts = num_parts;
    partition.element_part.resize(mesh.wedgeCount());
    partition.vertex_part.resize(mesh.vertexCount(), 0);

    if (mesh.wedgeCount() == 0) {
        return partition;
    }

    if (num_parts == 1) {
        std::fill(partition.element_part.begin(), partition.element_part.end(), 0);
        std::fill(partition.vertex_part.begin(), partition.vertex_part.end(), 0);
        return partition;
    }

    std::vector<idx_t> eptr(mesh.wedgeCount() + 1);
    std::vector<idx_t> eind(mesh.wedgeCount() * 6);
    for (std::size_t element = 0; element < mesh.wedgeCount(); ++element) {
        eptr[element] = static_cast<idx_t>(6 * element);
        for (int corner = 0; corner < 6; ++corner) {
            eind[6 * element + corner] = static_cast<idx_t>(mesh.wedges[element][corner]);
        }
    }
    eptr[mesh.wedgeCount()] = static_cast<idx_t>(eind.size());

    idx_t ne = static_cast<idx_t>(mesh.wedgeCount());
    idx_t nn = static_cast<idx_t>(mesh.vertexCount());
    idx_t ncommon = 4;
    idx_t nparts = static_cast<idx_t>(num_parts);
    idx_t objval = 0;

    std::vector<idx_t> epart(mesh.wedgeCount());
    std::vector<idx_t> npart(mesh.vertexCount());
    std::vector<idx_t> options(METIS_NOPTIONS);
    METIS_SetDefaultOptions(options.data());
    options[METIS_OPTION_NUMBERING] = 0;

    const int status = METIS_PartMeshNodal(
        &ne,
        &nn,
        eptr.data(),
        eind.data(),
        nullptr,
        nullptr,
        &nparts,
        nullptr,
        options.data(),
        &objval,
        epart.data(),
        npart.data());

    if (status != METIS_OK) {
        throw std::runtime_error("METIS_PartMeshNodal failed for wedge mesh");
    }

    for (std::size_t i = 0; i < epart.size(); ++i) {
        partition.element_part[i] = static_cast<int>(epart[i]);
    }
    for (std::size_t i = 0; i < npart.size(); ++i) {
        partition.vertex_part[i] = static_cast<int>(npart[i]);
    }

    return partition;
}

WedgeMesh3D MeshPartitioner::reorderByPartition(WedgeMesh3D mesh, MeshPartition& partition) {
    const int vertex_count = static_cast<int>(mesh.vertexCount());
    std::vector<int> order(vertex_count);
    for (int i = 0; i < vertex_count; ++i) {
        order[i] = i;
    }

    std::stable_sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        if (partition.vertex_part[lhs] != partition.vertex_part[rhs]) {
            return partition.vertex_part[lhs] < partition.vertex_part[rhs];
        }
        return lhs < rhs;
    });

    std::vector<int> old_to_new(vertex_count);
    std::vector<int> new_part(vertex_count);
    for (int new_index = 0; new_index < vertex_count; ++new_index) {
        const int old_index = order[new_index];
        old_to_new[old_index] = new_index;
        new_part[new_index] = partition.vertex_part[old_index];
    }

    WedgeMesh3D reordered;
    reordered.conductivity = mesh.conductivity;
    reordered.vertices.resize(mesh.vertices.size());
    for (int old_index = 0; old_index < vertex_count; ++old_index) {
        reordered.vertices[old_to_new[old_index]] = mesh.vertices[old_index];
    }

    reordered.wedges.reserve(mesh.wedges.size());
    for (const auto& wedge : mesh.wedges) {
        reordered.wedges.push_back({
            old_to_new[wedge[0]],
            old_to_new[wedge[1]],
            old_to_new[wedge[2]],
            old_to_new[wedge[3]],
            old_to_new[wedge[4]],
            old_to_new[wedge[5]],
        });
    }

    partition.vertex_part = std::move(new_part);
    partition.element_part.resize(mesh.wedgeCount());
    for (std::size_t element = 0; element < mesh.wedgeCount(); ++element) {
        const auto& wedge = mesh.wedges[element];
        partition.element_part[element] = partition.vertex_part[old_to_new[wedge[0]]];
    }

    return reordered;
}

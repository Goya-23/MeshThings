#ifndef _HEAT_FEM_H_
#define _HEAT_FEM_H_

#include <array>
#include <cstddef>
#include <map>
#include <vector>

#include <mpi.h>

#include "meshPartitioner.h"
#include "wedgeMesh.h"

struct LocalLinearSystem {
    int global_size = 0;
    int row_start = 0;
    int row_end = 0;
    int assembled_wedges = 0;
    std::vector<int> owned_vertices;
    std::vector<int> local_to_global;
    std::map<int, int> global_to_local;
    std::vector<int> row_indices;
    std::vector<int> row_ptr;
    std::vector<int> col_indices;
    std::vector<double> values;
    std::vector<double> rhs;
    std::vector<double> solution;
    std::vector<char> is_dirichlet;
};

class HeatFEMAssembler {
public:
    static LocalLinearSystem assemble(
        MPI_Comm comm,
        const WedgeMesh3D& mesh,
        const MeshPartition& partition);

    static bool isBoundaryVertex(const WedgeMesh3D& mesh, int vertex_id);

private:
    static double exactSolution(double x, double y, double z);
    static double sourceTerm(double x, double y, double z);
    static void addWedgeContribution(
        const WedgeMesh3D& mesh,
        std::size_t wedge_id,
        const std::array<int, 6>& wedge,
        std::map<int, std::map<int, double>>& matrix,
        std::map<int, double>& rhs);
    static void addTetrahedronContribution(
        const WedgeMesh3D& mesh,
        const std::array<int, 4>& tet,
        double conductivity,
        std::map<int, std::map<int, double>>& matrix,
        std::map<int, double>& rhs);
};

#endif

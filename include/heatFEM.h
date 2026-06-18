#ifndef _HEAT_FEM_H_
#define _HEAT_FEM_H_

#include <map>
#include <vector>

#include <mpi.h>

#include "meshPartitioner.h"
#include "triangulation2d.h"

struct LocalLinearSystem {
    int global_size = 0;
    int row_start = 0;
    int row_end = 0;
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
      const Triangulation2D& mesh,
      const MeshPartition& partition,
      double conductivity = 1.0);

private:
    static double exactSolution(double x, double y);
    static double sourceTerm(double x, double y);
    static bool isBoundaryVertex(const Triangulation2D& mesh, int vertex_id);
    static void addTriangleContribution(
        const Triangulation2D& mesh,
        const std::array<int, 3>& tri,
        double conductivity,
        std::map<int, std::map<int, double>>& matrix,
        std::map<int, double>& rhs);
};

#endif

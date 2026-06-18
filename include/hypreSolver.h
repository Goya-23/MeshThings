#ifndef _HYPRE_SOLVER_H_
#define _HYPRE_SOLVER_H_

#include <mpi.h>

#include "heatFEM.h"

class HypreSolver {
public:
    static std::vector<double> solve(MPI_Comm comm, const LocalLinearSystem& system);
    static std::vector<double> gatherSolution(MPI_Comm comm, const LocalLinearSystem& system,
                                              const std::vector<double>& local_solution);
};

#endif

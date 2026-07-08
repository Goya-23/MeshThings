#include "hypreSolver.h"

#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_parcsr_mv.h>
#include <HYPRE_utilities.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

namespace {

LocalLinearSystem gatherSystem(MPI_Comm comm, const LocalLinearSystem& local) {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    LocalLinearSystem global;
    global.global_size = local.global_size;

    if (rank == 0) {
        global.owned_vertices.resize(global.global_size);
        for (int i = 0; i < global.global_size; ++i) {
            global.owned_vertices[i] = i;
        }
        global.row_start = 0;
        global.row_end = global.global_size - 1;
        global.row_indices = global.owned_vertices;
        global.row_ptr.assign(global.global_size + 1, 0);
        global.rhs.assign(global.global_size, 0.0);
        global.solution.assign(global.global_size, 0.0);
    }

    const int local_rows = static_cast<int>(local.owned_vertices.size());
    const int nnz = static_cast<int>(local.col_indices.size());

    std::vector<int> row_counts(size);
    MPI_Gather(&local_rows, 1, MPI_INT, rank == 0 ? row_counts.data() : nullptr, 1, MPI_INT, 0, comm);

    std::vector<int> nnz_per_rank;
    if (rank == 0) {
        nnz_per_rank.resize(size);
    }
    MPI_Gather(&nnz, 1, MPI_INT, rank == 0 ? nnz_per_rank.data() : nullptr, 1, MPI_INT, 0, comm);

    if (rank != 0) {
        if (!local.owned_vertices.empty()) {
            MPI_Send(local.owned_vertices.data(), static_cast<int>(local.owned_vertices.size()), MPI_INT, 0, 0, comm);
            MPI_Send(local.row_ptr.data(), static_cast<int>(local.row_ptr.size()), MPI_INT, 0, 1, comm);
            MPI_Send(local.col_indices.data(), nnz, MPI_INT, 0, 2, comm);
            MPI_Send(const_cast<double*>(local.values.data()), nnz, MPI_DOUBLE, 0, 3, comm);
            MPI_Send(const_cast<double*>(local.rhs.data()), static_cast<int>(local.rhs.size()), MPI_DOUBLE, 0, 4, comm);
        }
        return global;
    }

    std::vector<std::vector<int>> all_rows(size);
    std::vector<std::vector<int>> all_row_ptr(size);
    std::vector<std::vector<int>> all_cols(size);
    std::vector<std::vector<double>> all_values(size);
    std::vector<std::vector<double>> all_rhs(size);

    for (int process = 0; process < size; ++process) {
        if (row_counts[process] == 0) {
            continue;
        }
        all_rows[process].resize(row_counts[process]);
        all_row_ptr[process].resize(row_counts[process] + 1);
        all_cols[process].resize(nnz_per_rank[process]);
        all_values[process].resize(nnz_per_rank[process]);
        all_rhs[process].resize(row_counts[process]);
        if (process == rank) {
            all_rows[process] = local.owned_vertices;
            all_row_ptr[process] = local.row_ptr;
            all_cols[process] = local.col_indices;
            all_values[process] = local.values;
            all_rhs[process] = local.rhs;
        } else {
            MPI_Recv(all_rows[process].data(), row_counts[process], MPI_INT, process, 0, comm, MPI_STATUS_IGNORE);
            MPI_Recv(all_row_ptr[process].data(), row_counts[process] + 1, MPI_INT, process, 1, comm, MPI_STATUS_IGNORE);
            MPI_Recv(all_cols[process].data(), nnz_per_rank[process], MPI_INT, process, 2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(all_values[process].data(), nnz_per_rank[process], MPI_DOUBLE, process, 3, comm, MPI_STATUS_IGNORE);
            MPI_Recv(all_rhs[process].data(), row_counts[process], MPI_DOUBLE, process, 4, comm, MPI_STATUS_IGNORE);
        }
    }

    std::vector<std::map<int, double>> matrix(global.global_size);
    for (int process = 0; process < size; ++process) {
        const int rows = row_counts[process];
        for (int row = 0; row < rows; ++row) {
            const int global_row = all_rows[process][row];
            global.rhs[global_row] = all_rhs[process][row];
            for (int entry = all_row_ptr[process][row]; entry < all_row_ptr[process][row + 1]; ++entry) {
                const int col = all_cols[process][entry];
                matrix[global_row][col] += all_values[process][entry];
            }
        }
    }

    global.row_ptr.assign(global.global_size + 1, 0);
    for (int row = 0; row < global.global_size; ++row) {
        global.row_ptr[row + 1] = global.row_ptr[row] + static_cast<int>(matrix[row].size());
    }

    global.col_indices.reserve(global.row_ptr.back());
    global.values.reserve(global.row_ptr.back());
    for (int row = 0; row < global.global_size; ++row) {
        for (const auto& entry : matrix[row]) {
            global.col_indices.push_back(entry.first);
            global.values.push_back(entry.second);
        }
    }

    return global;
}

std::vector<double> solveOnRank0(MPI_Comm comm, const LocalLinearSystem& system) {
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    const int local_rows = static_cast<int>(system.owned_vertices.size());
    std::vector<double> solution(local_rows, 0.0);
    if (rank != 0 || local_rows == 0) {
        return solution;
    }

    HYPRE_IJMatrix ij_matrix;
    HYPRE_ParCSRMatrix par_matrix;
    HYPRE_IJVector ij_rhs;
    HYPRE_IJVector ij_solution;
    HYPRE_ParVector par_rhs;
    HYPRE_ParVector par_solution;

    const HYPRE_BigInt row_start = system.row_start;
    const HYPRE_BigInt row_end = system.row_end;

    HYPRE_IJMatrixCreate(MPI_COMM_SELF, row_start, row_end, row_start, row_end, &ij_matrix);
    HYPRE_IJMatrixSetObjectType(ij_matrix, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(ij_matrix);

    for (int row = 0; row < local_rows; ++row) {
        const HYPRE_BigInt global_row = system.owned_vertices[row];
        HYPRE_Int num_cols = system.row_ptr[row + 1] - system.row_ptr[row];
        std::vector<HYPRE_BigInt> cols(num_cols);
        std::vector<double> values(num_cols);
        for (HYPRE_Int col = 0; col < num_cols; ++col) {
            cols[col] = system.col_indices[system.row_ptr[row] + col];
            values[col] = system.values[system.row_ptr[row] + col];
        }
        HYPRE_IJMatrixSetValues(ij_matrix, 1, &num_cols, &global_row, cols.data(), values.data());
    }

    HYPRE_IJMatrixAssemble(ij_matrix);
    HYPRE_IJMatrixGetObject(ij_matrix, reinterpret_cast<void**>(&par_matrix));

    HYPRE_IJVectorCreate(MPI_COMM_SELF, row_start, row_end, &ij_rhs);
    HYPRE_IJVectorSetObjectType(ij_rhs, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(ij_rhs);

    HYPRE_IJVectorCreate(MPI_COMM_SELF, row_start, row_end, &ij_solution);
    HYPRE_IJVectorSetObjectType(ij_solution, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(ij_solution);

    for (int row = 0; row < local_rows; ++row) {
        const HYPRE_BigInt global_row = system.owned_vertices[row];
        const double rhs_value = system.rhs[row];
        HYPRE_IJVectorSetValues(ij_rhs, 1, &global_row, &rhs_value);
        const double initial_guess = 0.0;
        HYPRE_IJVectorSetValues(ij_solution, 1, &global_row, &initial_guess);
    }

    HYPRE_IJVectorAssemble(ij_rhs);
    HYPRE_IJVectorAssemble(ij_solution);
    HYPRE_IJVectorGetObject(ij_rhs, reinterpret_cast<void**>(&par_rhs));
    HYPRE_IJVectorGetObject(ij_solution, reinterpret_cast<void**>(&par_solution));

    HYPRE_Solver solver;
    HYPRE_Solver precond;
    HYPRE_ParCSRPCGCreate(MPI_COMM_SELF, &solver);
    HYPRE_ParCSRPCGSetMaxIter(solver, 2000);
    HYPRE_ParCSRPCGSetTol(solver, 1.0e-10);
    HYPRE_ParCSRPCGSetTwoNorm(solver, 1);
    HYPRE_ParCSRPCGSetPrintLevel(solver, 0);
    HYPRE_BoomerAMGCreate(&precond);
    HYPRE_BoomerAMGSetPrintLevel(precond, 0);
    HYPRE_BoomerAMGSetTol(precond, 0.0);
    HYPRE_BoomerAMGSetMaxIter(precond, 1);
    HYPRE_ParCSRPCGSetPrecond(
        solver,
        HYPRE_BoomerAMGSolve,
        HYPRE_BoomerAMGSetup,
        precond);
    HYPRE_ParCSRPCGSetup(solver, par_matrix, par_rhs, par_solution);
    HYPRE_ParCSRPCGSolve(solver, par_matrix, par_rhs, par_solution);

    HYPRE_Int iterations = 0;
    HYPRE_Real residual = 0.0;
    HYPRE_ParCSRPCGGetNumIterations(solver, &iterations);
    HYPRE_ParCSRPCGGetFinalRelativeResidualNorm(solver, &residual);
    std::cout << "HYPRE PCG iterations: " << iterations << ", residual: " << residual << "\n";

    for (int row = 0; row < local_rows; ++row) {
        HYPRE_BigInt global_row = system.owned_vertices[row];
        double value = 0.0;
        HYPRE_ParVectorGetValues(par_solution, 1, &global_row, &value);
        solution[row] = value;
    }

    HYPRE_ParCSRPCGDestroy(solver);
    HYPRE_BoomerAMGDestroy(precond);
    HYPRE_IJMatrixDestroy(ij_matrix);
    HYPRE_IJVectorDestroy(ij_rhs);
    HYPRE_IJVectorDestroy(ij_solution);
    return solution;
}

}  // namespace

std::vector<double> HypreSolver::solve(MPI_Comm comm, const LocalLinearSystem& system) {
    LocalLinearSystem global_system = gatherSystem(comm, system);
    std::vector<double> global_solution = solveOnRank0(comm, global_system);

    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    const int global_size = system.global_size;
    std::vector<double> broadcast_solution(global_size, 0.0);
    if (rank == 0) {
        broadcast_solution = global_solution;
    }
    MPI_Bcast(broadcast_solution.data(), global_size, MPI_DOUBLE, 0, comm);

    std::vector<double> local_solution;
    local_solution.reserve(system.owned_vertices.size());
    for (int vertex : system.owned_vertices) {
        local_solution.push_back(broadcast_solution[vertex]);
    }
    return local_solution;
}

std::vector<double> HypreSolver::gatherSolution(
    MPI_Comm comm, const LocalLinearSystem& system, const std::vector<double>& local_solution) {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::vector<double> global_solution(system.global_size, 0.0);
    for (std::size_t row = 0; row < system.owned_vertices.size(); ++row) {
        global_solution[system.owned_vertices[row]] = local_solution[row];
    }

    std::vector<double> reduced(system.global_size, 0.0);
    MPI_Allreduce(global_solution.data(), reduced.data(), system.global_size, MPI_DOUBLE, MPI_SUM, comm);
    return reduced;
}

#pragma once

#include <string>
#include <memory>
#include <vector>

#include <Eigen/Core>
#include <ortools/linear_solver/linear_solver.h>

namespace BRY {

namespace ort = operations_research;

class LPSolver {
    public:
        struct Result {
            /// @brief Optimizer result
            ort::MPSolver::ResultStatus status;

            /// @brief Probability of safety
            bry_float_t p_safe;
        
            /// @brief eta (init set constraint) and gamma (expected increase constraint)
            bry_float_t eta, gamma;

            /// @brief Optimal barrier coefficients
            Vector b_values;

            double comp_time;
        };

    public:
        /// @brief Create a LP solver using OR-Tools
        /// @param solver_id Solver ID passed to OR-Tools MPSolver
        LPSolver(const std::string& solver_id);

        void setConstraintMatrices(const Matrix& A, const Vector& b);

        /// @brief Give a hint to the optimizer using the trivial barrier (B(x) = 1, eta = 1, gamma = 0)
        void setTrivialBarrierHint();

        /// @brief Set a time limit for the solver
        /// @param time_limit_ms Time limit in milliseconds
        void setTimeLimit(int64_t time_limit_ms);

        /// @brief Solve the optimization problem for a given time horizon
        /// @param time_horizon Number of timesteps to verify 
        Result solve(uint32_t time_horizon);

        /// @brief Get the solution vector `v` that solves `Av >= b` once the LP has been solved
        /// @return Solution vector in the form (b_0, ..., b_m, eta, gamma)
        Vector getSolnVector() const;
    private:
        BRY_INL bry_int_t nMonoms() const {return m_b.size();}

    private:
        std::unique_ptr<ort::MPSolver> m_solver;
        std::vector<ort::MPVariable*> m_b;
        ort::MPVariable* m_eta;
        ort::MPVariable* m_gamma;
        ort::MPObjective* m_objective = nullptr;
        ort::MPSolver::ResultStatus m_result_status;

        bool m_constraints_set = false;
        bry_float_t m_inf;
};

}

#include "impl/LPSolver_impl.hpp"
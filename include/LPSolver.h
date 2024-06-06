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
        };

    public:
        /// @brief Create a LP solver using OR-Tools
        /// @param solver_id Solver ID passed to OR-Tools MPSolver
        /// @param n_monoms Number of monomials in barrier
        LPSolver(const std::string& solver_id, bry_int_t n_monoms);

        /// @brief Set the constraint for non-negativity over the workspace (B(x) >= 0)
        /// @param b_coeffs Coefficient matrix for b vector
        void addWorkspaceConstraint(const Matrix& b_coeffs);

        /// @brief Add constraint for an intial set region (B(x) - eta <= 0)
        /// @param b_coeffs Coefficient matrix for b vector
        /// @param eta_coeffs Coefficient vector for eta (number of elements must match rows in `b_coeffs`)
        void addInitialSetConstraint(const Matrix& b_coeffs, const Vector& eta_coeffs);

        /// @brief Add constraint for an unsafe set region (B(x) > 1)
        /// @param b_coeffs Coefficient matrix for b vector
        /// @param lower_bound Vector of lower bound scalars (RHS)
        void addUnsafeSetConstraint(const Matrix& b_coeffs, const Vector& lower_bound);

        /// @brief Add constraint for a safe set region (E[B(f(x) + v)] - B(x) - gamma <= 0)
        /// @param b_coeffs Coefficient matrix for b vector
        /// @param gamma_coeffs Coefficient vector for gamma (number of elements must match rows in `b_coeffs`)
        void addSafeSetConstraint(const Matrix& b_coeffs, const Vector& gamma_coeffs);

        /// @brief Give a hint to the optimizer using the trivial barrier (B(x) = 1, eta = 1, gamma = 0)
        void setTrivialBarrierHint();

        /// @brief Set a time limit for the solver
        /// @param time_limit_ms Time limit in milliseconds
        void setTimeLimit(int64_t time_limit_ms);

        /* TODO */
        //void enforceBoundaryConstraints(const Matrix& constraint);

        /// @brief Solve the optimization problem for a given time horizon
        /// @param time_horizon Number of timesteps to verify 
        Result solve(uint32_t time_horizon);

    private:
        std::unique_ptr<ort::MPSolver> m_solver;
        bry_int_t m_n_monoms;
        bry_float_t m_inf;
        std::vector<ort::MPVariable*> m_b;
        ort::MPVariable* m_eta;
        ort::MPVariable* m_gamma;
};

}

#include "impl/LPSolver_impl.hpp"
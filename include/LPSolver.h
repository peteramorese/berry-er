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
        
            /// @brief Optimal barrier coefficients
            Eigen::VectorXd beta_values;
        };

    public:
        /// @brief Create a LP solver using OR-Tools
        /// @param solver_id Solver ID passed to OR-Tools MPSolver
        /// @param n_monoms Number of monomials in barrier
        LPSolver(const std::string& solver_id, bry_deg_t n_monoms);

        /// @brief Add constraint for an intial set region (B(x) - eta <= 0)
        /// @param beta_coeffs Coefficient matrix for beta vector
        /// @param eta_coeffs Coefficient vector for eta (number of elements must match rows in `beta_coeffs`)
        void addInitialSetConstraint(const Eigen::MatrixXd& beta_coeffs, const Eigen::VectorXd& eta_coeffs);

        /// @brief Add constraint for an unsafe set region (B(x) > 1)
        /// @param beta_coeffs Coefficient matrix for beta vector
        /// @param lower_bound Vector of lower bound scalars (RHS)
        void addUnsafeSetConstraint(const Eigen::MatrixXd& beta_coeffs, const Eigen::VectorXd& lower_bound);

        /// @brief Add constraint for a safe set region (E[B(f(x) + v)] - B(x) - gamma <= 0)
        /// @param beta_coeffs Coefficient matrix for beta vector
        /// @param gamma_coeffs Coefficient vector for gamma (number of elements must match rows in `beta_coeffs`)
        void addSafeSetConstraint(const Eigen::MatrixXd& beta_coeffs, const Eigen::VectorXd& gamma_coeffs);

        /* TODO */
        //void enforceBoundaryConstraints(const Eigen::MatrixXd& constraint);

        /// @brief Solve the optimization problem for a given time horizon
        /// @param time_horizon Number of timesteps to verify 
        Result solve(uint32_t time_horizon);

    private:
        std::unique_ptr<ort::MPSolver> m_solver;
        bry_deg_t m_n_monoms;
        bry_float_t m_inf;
        std::vector<ort::MPVariable*> m_beta;
        ort::MPVariable* m_eta;
        ort::MPVariable* m_gamma;
};

}

#include "impl/LPSolver_impl.hpp"
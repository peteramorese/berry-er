#pragma once

#include "LPSolver.h"

#include "berry/Options.h"

#include <ortools/linear_solver/linear_solver.h>

BRY::LPSolver::LPSolver(const std::string& solver_id, bry_deg_t n_monoms)
    : m_solver(MPSolver::CreateSolver("SCIP"))
    , m_n_variables(n_monoms + 2)
{
    if (!m_solver) {
        LOG(FATAL) << "Solver unavailable.";
        ASSERT(false, solver_id << " is not available");
    } else {
        INFO("Solver '" << solver_id << "' found!")
        INFO("Constructing LP with " << n_variables << " variables");
        m_inf = m_solver->infinity();
        m_solver->MakeNumVarArray(n_monoms, 0.0, m_inf, "beta", m_beta);
        m_eta = m_solver->MakeNumVar(0.0, m_inf, "eta");
        m_gamma = m_solver->MakeNumVar(0.0, m_inf, "gamma");
    }
}


void BRY::LPSolver::addInitialSetConstraint(const Eigen::MatrixXd& beta_coeffs, const Eigen::VectorXd& eta_coeffs) {
    bry_idx_t rows = eta_coeffs.size();
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT(beta_coeffs.rows() == rows, "Number of rows in `beta_coeffs` does not match elements in `eta_coeffs`");
    ASSERT(beta_coeffs.cols() == m_n_monoms, "Number of columns does not match number of beta monomials");
#endif
    for (bry_idx_t i = 0; i < rows; ++i) {
        ort::MPConstraint* row_constraint = m_solver->MakeRowConstraint(0.0, m_inf, "");
        for (bry_idx_t j = 0; j < constraint_coeffs.cols(); ++j) {
            row_constraint->SetCoefficient(m_beta[j], constraint(i, j));
        }
        row_constraint->SetCoefficient(m_eta, eta_coeffs(i));
        // Gamma does not appear in this constraint
        row_constraint->SetCoefficient(m_gamma, 0.0);
    }
    INFO("[Initial set] Added " << rows << " constraints");
}

void BRY::LPSolver::addUnsafeSetConstraint(const Eigen::MatrixXd& beta_coeffs, const Eigen::VectorXd& lower_bound) {
    bry_idx_t rows = lower_bound.size();
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT(beta_coeffs.rows() == rows, "Number of rows in `beta_coeffs` does not match elements in `lower_bound`");
    ASSERT(beta_coeffs.cols() == m_n_monoms, "Number of columns does not match number of beta monomials");
#endif
    for (bry_idx_t i = 0; i < rows; ++i) {
        ort::MPConstraint* row_constraint = m_solver->MakeRowConstraint(lower_bound(i), m_inf, "");
        for (bry_idx_t j = 0; j < constraint_coeffs.cols(); ++j) {
            row_constraint->SetCoefficient(m_beta[j], constraint(i, j));
        }
        // Neither eta nor gamma appear in this constraint
        row_constraint->SetCoefficient(m_eta, 0.0);
        row_constraint->SetCoefficient(m_gamma, 0.0);
    }
    INFO("[Unsafe set] Added " << rows << " constraints");
}

void BRY::LPSolver::addSafeSetConstraint(const Eigen::MatrixXd& beta_coeffs, const Eigen::VectorXd& gamma_coeffs) {
    bry_idx_t rows = gamma_coeffs.size();
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT(beta_coeffs.rows() == rows, "Number of rows in `beta_coeffs` does not match elements in `gamma_coeffs`");
    ASSERT(beta_coeffs.cols() == m_n_monoms, "Number of columns does not match number of beta monomials");
#endif
    for (bry_idx_t i = 0; i < rows; ++i) {
        ort::MPConstraint* row_constraint = m_solver->MakeRowConstraint(lower_bound(i), m_inf, "");
        for (bry_idx_t j = 0; j < constraint_coeffs.cols(); ++j) {
            row_constraint->SetCoefficient(m_beta[j], constraint(i, j));
        }
        // Eta does not appear in this constraint
        row_constraint->SetCoefficient(m_eta, 0.0);
        row_constraint->SetCoefficient(m_gamma, gamma_coeffs(i));
    }
    INFO("[Safe set] Added " << rows << " constraints");
}

BRY::LPSolver::Result BRY::LPSolver::solve(uint32_t time_horizon) {
    ort::MPObjective* objective = m_solver->MutableObjective();

    // Beta do not contribute to objective
    for (ort::MPVariable* beta_i : m_beta)
        objective->SetCoefficient(beta_i, 0.0);

    objective->SetCoefficient(m_eta, 1.0);
    objective->SetCoefficient(m_gamma, static_cast<bry_float_t>(time_horizon));
    objective->SetMinimization()

    const MPSolver::ResultStatus result_status = m_solver->Solve();

    Eigen::VectorXd beta_values(m_n_monoms);
    for (bry_idx_t i = 0; i < m_n_monoms; ++i) {
        beta_values(i) = m_beta[i]->solution_value();
    }

    return Result{result_status, 1.0 - objective->Value(), std::move(beta_values)};
}
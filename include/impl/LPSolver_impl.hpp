#pragma once

#include "LPSolver.h"

#include "Time.h"

#include "berry/Options.h"

#include <ortools/linear_solver/linear_solver.pb.h>

BRY::LPSolver::LPSolver(const std::string& solver_id)
    : m_solver(ort::MPSolver::CreateSolver(solver_id))
    , m_constraints_set(false)
{
    if (!m_solver) {
        LOG(FATAL) << "Solver unavailable.";
        ASSERT(false, solver_id << " is not available");
    } else {
        INFO("Solver '" << solver_id << "' found!");
        m_inf = m_solver->infinity();
    }
}

void BRY::LPSolver::setConstraintMatrices(const Matrix& A, const Vector& b) {
    if (m_constraints_set) {
        ERROR("Constraints matrices have already been set");
        return;
    }
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT(A.rows() == b.size(), "Number of rows in `A` does not match size of `b`");
#endif

    INFO("Constructing LP with " << A.cols() << " variables");
    m_solver->MakeNumVarArray(A.cols() - 2, -m_inf, m_inf, "b", &m_b);

#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT(A.cols() == nMonoms() + 2, "Number of columns does not match number of variables");
#endif

    m_eta = m_solver->MakeNumVar(0.0, m_inf, "eta");
    m_gamma = m_solver->MakeNumVar(0.0, m_inf, "gamma");

    for (bry_int_t i = 0; i < A.rows(); ++i) {
        ort::MPConstraint* row_constraint = m_solver->MakeRowConstraint(b(i), m_inf, "");
        for (bry_int_t j = 0; j < A.cols() - 2; ++j) {
            row_constraint->SetCoefficient(m_b[j], A(i, j));
        }
        row_constraint->SetCoefficient(m_eta, A(i, A.cols() - 2));
        row_constraint->SetCoefficient(m_gamma, A(i, A.cols() - 1));
    }
    INFO("Added " << A.rows() << " constraints");
}


void BRY::LPSolver::setTrivialBarrierHint() {
    std::vector<std::pair<const ort::MPVariable*, double>> hint;
    hint.reserve(nMonoms() + 2);
    for (const ort::MPVariable* b_i : m_b) {
        hint.push_back({b_i, 0.0});
    }
    hint.front() = {m_b[0], 1.0};
    hint.push_back({m_eta, 1.0});
    hint.push_back({m_gamma, 0.0});
    m_solver->SetHint(hint);
}

void BRY::LPSolver::setTimeLimit(int64_t time_limit_ms) {
    m_solver->set_time_limit(time_limit_ms);
}

BRY::LPSolver::Result BRY::LPSolver::solve(uint32_t time_horizon) {
    ort::MPObjective* objective = m_solver->MutableObjective();

    // Beta do not contribute to objective
    for (ort::MPVariable* b_i : m_b)
        objective->SetCoefficient(b_i, 0.0);

    objective->SetCoefficient(m_eta, 1.0);
    objective->SetCoefficient(m_gamma, static_cast<bry_float_t>(time_horizon));
    objective->SetMinimization();

    Timer t("synthesis");
    const ort::MPSolver::ResultStatus result_status = m_solver->Solve();

    Vector beta_values(nMonoms());
    for (bry_int_t i = 0; i < nMonoms(); ++i) {
        beta_values(i) = m_b[i]->solution_value();
    }

    return Result{result_status, 1.0 - objective->Value(), m_eta->solution_value(), m_gamma->solution_value(), std::move(beta_values), t.now(TimeUnit::s)};
}
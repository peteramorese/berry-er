#pragma once

#include "Synthesis.h"

#include "berry/BernsteinTransform.h"

/* PolyDynamicsSynthesizer */

template <std::size_t DIM>
BRY::PolyDynamicsSynthesizer<DIM>::PolyDynamicsSynthesizer(const std::shared_ptr<PolynomialDynamics<DIM>>& dynamics, const std::shared_ptr<Additive2ndMomentNoise<DIM>>& noise, bry_deg_t barrier_deg, const std::string& solver_id)
    : m_dynamics(dynamics)
    , m_solver(new LPSolver(solver_id, pow(barrier_deg + 1, DIM)))
    , m_barrier_deg(barrier_deg)
    , m_noise(noise)
{}

template <std::size_t DIM>
void BRY::PolyDynamicsSynthesizer<DIM>::setProblem(const std::shared_ptr<SynthesisProblem<DIM>>& problem) {
    m_problem = problem;
}

template <std::size_t DIM>
std::pair<Eigen::MatrixXd, Eigen::VectorXd> BRY::PolyDynamicsSynthesizer<DIM>::getConstraintMatrices(bry_deg_t degree_increase) const {
    Eigen::MatrixXd Phi_m = BernsteinBasisTransform<DIM>::pwrToBernMatrix(m_barrier_deg, degree_increase);
    Eigen::MatrixXd Phi_inv_m = BernsteinBasisTransform<DIM>::bernToPwrMatrix(m_barrier_deg);

    // Workspace
    bry_deg_t n_cols = Phi_m.cols() + 2;
    Eigen::MatrixXd ws_beta_coeffs = Phi_m * m_problem->workspace.transformationMatrix(m_barrier_deg) * Phi_inv_m;
    Eigen::MatrixXd ws_coeffs(ws_beta_coeffs.rows(), n_cols);
    //           beta            eta                                           gamma
    ws_coeffs << ws_beta_coeffs, Eigen::VectorXd::Zero(ws_beta_coeffs.rows()), Eigen::VectorXd::Zero(ws_beta_coeffs.rows());
    Eigen::VectorXd ws_lower_bound = Eigen::VectorXd::Zero(ws_beta_coeffs.rows());

    // Initial sets
    if (m_problem->init_sets.empty())
        WARN("No initial sets were provided");
    Eigen::MatrixXd init_coeffs(m_problem->init_sets.size() * Phi_m.rows(), n_cols);
    Eigen::VectorXd init_lower_bound = Eigen::VectorXd::Zero(init_coeffs.rows());
    bry_idx_t i = 0;
    for (const HyperRectangle<DIM>& set : m_problem->init_sets) {
        Eigen::MatrixXd beta_coeffs = -Phi_m * set.transformationMatrix(m_barrier_deg) * Phi_inv_m;
        //DEBUG("init beta coeffs \n" << beta_coeffs);
        Eigen::MatrixXd coeffs(beta_coeffs.rows(), n_cols);
        //        beta         eta                                        gamma
        coeffs << beta_coeffs, Eigen::VectorXd::Ones(beta_coeffs.rows()), Eigen::VectorXd::Zero(beta_coeffs.rows());
        init_coeffs.block(Phi_m.rows() * i++, 0, Phi_m.rows(), n_cols) = coeffs;
    }

    // Unsafe sets
    if (m_problem->unsafe_sets.empty())
        WARN("No unsafe sets were provided");
    Eigen::MatrixXd unsafe_coeffs(m_problem->unsafe_sets.size() * Phi_m.rows(), n_cols);
    Eigen::VectorXd unsafe_lower_bound = Eigen::VectorXd::Ones(unsafe_coeffs.rows());
    i = 0;
    for (const HyperRectangle<DIM>& set : m_problem->unsafe_sets) {
        Eigen::MatrixXd beta_coeffs = Phi_m * set.transformationMatrix(m_barrier_deg) * Phi_inv_m;
        //DEBUG("unsafe beta coeffs \n" << beta_coeffs);
        Eigen::MatrixXd coeffs(beta_coeffs.rows(), n_cols);
        //        beta         eta                                        gamma
        coeffs << beta_coeffs, Eigen::VectorXd::Zero(beta_coeffs.rows()), Eigen::VectorXd::Zero(beta_coeffs.rows());
        unsafe_coeffs.block(Phi_m.rows() * i++, 0, Phi_m.rows(), n_cols) = coeffs;
    }

    // Safe sets
    if (m_problem->safe_sets.empty())
        WARN("No safe sets were provided");
    Eigen::MatrixXd F_expec_Gamma = m_dynamics->dynamicsPowerMatrix(m_barrier_deg) * m_noise->additiveNoiseMatrix(m_barrier_deg);
    bry_deg_t p = m_dynamics->composedDegree(m_barrier_deg);
    Eigen::MatrixXd Phi_p = BernsteinBasisTransform<DIM>::pwrToBernMatrix(p, degree_increase);
    Eigen::VectorXd gamma_coeffs = Eigen::VectorXd::Ones(Phi_p.cols());
    ASSERT(F_expec_Gamma.rows() == Phi_p.cols(), "Dimension mismatch between F and Phi (p)");

    Eigen::MatrixXd safe_coeffs(m_problem->safe_sets.size() * Phi_p.rows(), n_cols);
    Eigen::VectorXd safe_lower_bound = Eigen::VectorXd::Zero(safe_coeffs.rows());
    i = 0;
    for (const HyperRectangle<DIM>& set : m_problem->safe_sets) {
        Eigen::MatrixXd beta_coeffs = 
            -Phi_p * set.transformationMatrix(p) * (F_expec_Gamma - Eigen::MatrixXd::Identity(F_expec_Gamma.rows(), F_expec_Gamma.cols())) * Phi_inv_m;
        
        //DEBUG("safe beta coeffs \n" << beta_coeffs);
        Eigen::MatrixXd coeffs(beta_coeffs.rows(), n_cols);
        //        beta         eta                                        gamma
        coeffs << beta_coeffs, Eigen::VectorXd::Zero(beta_coeffs.rows()), Eigen::VectorXd::Ones(beta_coeffs.rows());
        safe_coeffs.block(Phi_p.rows() * i++, 0, Phi_p.rows(), n_cols) = coeffs;
    }

    Eigen::MatrixXd A(ws_coeffs.rows() + init_coeffs.rows() + unsafe_coeffs.rows() + safe_coeffs.rows(), n_cols);
    A << ws_coeffs, init_coeffs, unsafe_coeffs, safe_coeffs;

    Eigen::VectorXd b(ws_lower_bound.size() + init_lower_bound.size() + unsafe_lower_bound.size() + safe_lower_bound.size());
    b << ws_lower_bound, init_lower_bound, unsafe_lower_bound, safe_lower_bound;
    return std::make_pair(std::move(A), std::move(b));
}

template <std::size_t DIM>
void BRY::PolyDynamicsSynthesizer<DIM>::initialize(bry_deg_t degree_increase) {
    Eigen::MatrixXd Phi_m = BernsteinBasisTransform<DIM>::pwrToBernMatrix(m_barrier_deg, degree_increase);
    Eigen::MatrixXd Phi_inv_m = BernsteinBasisTransform<DIM>::bernToPwrMatrix(m_barrier_deg);
    //DEBUG("Phi inv m: \n" << Phi_inv_m);
    //DEBUG("Phi m: \n" << Phi_m);

    // Workspace
    {
        Eigen::MatrixXd beta_coeffs = Phi_m * m_problem->workspace.transformationMatrix(m_barrier_deg) * Phi_inv_m;
        m_solver->setWorkspaceConstraint(beta_coeffs);
    }

    // Initial sets
    if (m_problem->init_sets.empty())
        WARN("No initial sets were provided");
    Eigen::VectorXd eta_coeffs = Eigen::VectorXd::Ones(Phi_m.rows());
    for (const HyperRectangle<DIM>& set : m_problem->init_sets) {
        Eigen::MatrixXd beta_coeffs = -Phi_m * set.transformationMatrix(m_barrier_deg) * Phi_inv_m;
        m_solver->addInitialSetConstraint(beta_coeffs, eta_coeffs);
    }

    // Unsafe sets
    if (m_problem->unsafe_sets.empty())
        WARN("No unsafe sets were provided");
    Eigen::VectorXd lower_bound = Eigen::VectorXd::Ones(Phi_m.rows());
    for (const HyperRectangle<DIM>& set : m_problem->unsafe_sets) {
        Eigen::MatrixXd beta_coeffs = Phi_m * set.transformationMatrix(m_barrier_deg) * Phi_inv_m;
        m_solver->addUnsafeSetConstraint(beta_coeffs, lower_bound);
    }

    // Safe sets
    if (m_problem->safe_sets.empty())
        WARN("No safe sets were provided");
    Eigen::MatrixXd F_expec_Gamma = m_dynamics->dynamicsPowerMatrix(m_barrier_deg) * m_noise->additiveNoiseMatrix(m_barrier_deg);
    //DEBUG("F \n" << m_dynamics->dynamicsPowerMatrix(m_barrier_deg));
    //DEBUG("F expec Gamma: \n" << F_expec_Gamma);
    bry_deg_t p = m_dynamics->composedDegree(m_barrier_deg);
    Eigen::MatrixXd Phi_p = BernsteinBasisTransform<DIM>::pwrToBernMatrix(p, degree_increase);
    //Eigen::VectorXd gamma_coeffs = Phi_inv_p * Eigen::VectorXd::Ones(Phi_inv_p.cols());
    Eigen::VectorXd gamma_coeffs = Eigen::VectorXd::Ones(Phi_p.rows());
    ASSERT(F_expec_Gamma.rows() == Phi_p.cols(), "Dimension mismatch between F and Phi (p)");

    for (const HyperRectangle<DIM>& set : m_problem->safe_sets) {
        //Eigen::MatrixXd beta_coeffs = 
        //    -Phi_p * (F_expec_Gamma - Eigen::MatrixXd::Identity(F_expec_Gamma.rows(), F_expec_Gamma.cols())) 
        //    * set.transformationMatrix(m_barrier_deg) * Phi_inv_m;

        Eigen::MatrixXd beta_coeffs = 
            -Phi_p * set.transformationMatrix(p) * (F_expec_Gamma - Eigen::MatrixXd::Identity(F_expec_Gamma.rows(), F_expec_Gamma.cols())) * Phi_inv_m;

        //DEBUG("Safe set beta coeffs: \n" << beta_coeffs);
        //DEBUG("transformation matrix: \n" << set.transformationMatrix(m_barrier_deg));
        //DEBUG("matrix: \n" << (F_expec_Gamma - Eigen::MatrixXd::Identity(F_expec_Gamma.rows(), F_expec_Gamma.cols())) * set.transformationMatrix(m_barrier_deg));
        m_solver->addSafeSetConstraint(beta_coeffs, gamma_coeffs);
    }
    
    m_solver->setTrivialBarrierHint();

    m_initialized = true;
}

template <std::size_t DIM>
void BRY::PolyDynamicsSynthesizer<DIM>::setTimeLimit(int64_t time_limit_ms) {
    m_solver->setTimeLimit(time_limit_ms);
}

template <std::size_t DIM>
BRY::SynthesisSolution<DIM> BRY::PolyDynamicsSynthesizer<DIM>::synthesize() {
    if (!m_initialized) {
        ERROR("Synthesizer not initialized (call initialize() before synthesize())");
        return SynthesisSolution<DIM>{};
    }
    LPSolver::Result solver_result = m_solver->solve(m_problem->time_horizon);

    SynthesisSolution<DIM> solution;
    solution.p_safe = solver_result.p_safe;
    solution.eta = solver_result.eta;
    solution.gamma = solver_result.gamma;

    //DEBUG("min beta: " << solver_result.beta_values.minCoeff());
    //DEBUG("beta values: " << solver_result.beta_values.transpose());

    // TODO
    solution.success = true;

    solution.certificate.reset(new Polynomial<DIM, Basis::Bernstein>(solver_result.beta_values));
    return solution;
}

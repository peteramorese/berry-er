#pragma once

#include "Synthesis.h"

#include "berry/BernsteinTransform.h"

/* Synthesizer */

template <std::size_t DIM>
void BRY::Synthesizer<DIM>::setInitialSets(const std::vector<HyperRectangle<DIM>>& sets) {
    m_init_sets = sets;
}

template <std::size_t DIM>
void BRY::Synthesizer<DIM>::setInitialSets(std::vector<HyperRectangle<DIM>>&& sets) {
    m_init_sets = std::move(sets);
}

template <std::size_t DIM>
void BRY::Synthesizer<DIM>::setSafeSets(const std::vector<HyperRectangle<DIM>>& sets) {
    m_safe_sets = sets;
}

template <std::size_t DIM>
void BRY::Synthesizer<DIM>::setSafeSets(std::vector<HyperRectangle<DIM>>&& sets) {
    m_safe_sets = std::move(sets);
}

template <std::size_t DIM>
void BRY::Synthesizer<DIM>::setUnsafeSets(const std::vector<HyperRectangle<DIM>>& sets) {
    m_unsafe_sets = sets;
}

template <std::size_t DIM>
void BRY::Synthesizer<DIM>::setUnsafeSets(std::vector<HyperRectangle<DIM>>&& sets) {
    m_unsafe_sets = std::move(sets);
}

template <std::size_t DIM>
void BRY::Synthesizer<DIM>::setWorkspace(const HyperRectangle<DIM>& set) {
    m_workspace = set;
}


/* PolyDynamicsSynthesizer */

template <std::size_t DIM>
BRY::PolyDynamicsSynthesizer<DIM>::PolyDynamicsSynthesizer(const std::shared_ptr<PolynomialDynamics<DIM>>& dynamics, const std::shared_ptr<Additive2ndMomentNoise<DIM>>& noise, bry_deg_t barrier_deg)
    : m_dynamics(dynamics)
    //, m_solver(new LPSolver("SCIP", pow(barrier_deg + 1, DIM)))
    , m_solver(new LPSolver("PDLP", pow(barrier_deg + 1, DIM)))
    , m_barrier_deg(barrier_deg)
    , m_noise(noise)
{}

template <std::size_t DIM>
void BRY::PolyDynamicsSynthesizer<DIM>::initialize() {
    Eigen::MatrixXd Phi_m = BernsteinBasisTransform<DIM>::pwrToBernMatrix(m_barrier_deg);
    Eigen::MatrixXd Phi_inv_m = BernsteinBasisTransform<DIM>::bernToPwrMatrix(m_barrier_deg);
    //DEBUG("Phi inv m: \n" << Phi_inv_m);
    //DEBUG("Phi m: \n" << Phi_m);

    // Workspace
    {
        Eigen::MatrixXd beta_coeffs = Phi_m * this->m_workspace.transformationMatrix(m_barrier_deg) * Phi_inv_m;
        m_solver->setWorkspaceConstraint(beta_coeffs);
    }

    // Initial sets
    if (this->m_init_sets.empty())
        WARN("No initial sets were provided");
    Eigen::VectorXd eta_coeffs = Eigen::VectorXd::Ones(Phi_m.cols());
    for (const HyperRectangle<DIM>& set : this->m_init_sets) {
        Eigen::MatrixXd beta_coeffs = -Phi_m * set.transformationMatrix(m_barrier_deg) * Phi_inv_m;
        m_solver->addInitialSetConstraint(beta_coeffs, eta_coeffs);
    }

    // Unsafe sets
    if (this->m_unsafe_sets.empty())
        WARN("No unsafe sets were provided");
    Eigen::VectorXd lower_bound = Eigen::VectorXd::Ones(Phi_m.cols());
    for (const HyperRectangle<DIM>& set : this->m_unsafe_sets) {
        Eigen::MatrixXd beta_coeffs = Phi_m * set.transformationMatrix(m_barrier_deg) * Phi_inv_m;
        m_solver->addUnsafeSetConstraint(beta_coeffs, eta_coeffs);
    }

    // Safe sets
    if (this->m_safe_sets.empty())
        WARN("No safe sets were provided");
    Eigen::MatrixXd F_expec_Gamma = m_dynamics->dynamicsPowerMatrix(m_barrier_deg) * m_noise->additiveNoiseMatrix(m_barrier_deg);
    //DEBUG("F \n" << m_dynamics->dynamicsPowerMatrix(m_barrier_deg));
    //DEBUG("F expec Gamma: \n" << F_expec_Gamma);
    Eigen::MatrixXd Phi_inv_p = BernsteinBasisTransform<DIM>::pwrToBernMatrix(m_dynamics->composedDegree(m_barrier_deg));
    //Eigen::VectorXd gamma_coeffs = Phi_inv_p * Eigen::VectorXd::Ones(Phi_inv_p.cols());
    Eigen::VectorXd gamma_coeffs = Eigen::VectorXd::Ones(Phi_inv_p.cols());
    ASSERT(F_expec_Gamma.rows() == Phi_inv_p.cols(), "Dimension mismatch between F and Phi inv (p)");

    for (const HyperRectangle<DIM>& set : this->m_safe_sets) {
        DEBUG("b4 comp");
        Eigen::MatrixXd beta_coeffs = 
            -Phi_inv_p * (F_expec_Gamma - Eigen::MatrixXd::Identity(F_expec_Gamma.rows(), F_expec_Gamma.cols())) 
            * set.transformationMatrix(m_barrier_deg) * Phi_inv_m;
        DEBUG("af comp");
        //DEBUG("Safe set beta coeffs: \n" << beta_coeffs);
        //DEBUG("transformation matrix: \n" << set.transformationMatrix(m_barrier_deg));
        //DEBUG("matrix: \n" << (F_expec_Gamma - Eigen::MatrixXd::Identity(F_expec_Gamma.rows(), F_expec_Gamma.cols())) * set.transformationMatrix(m_barrier_deg));
        DEBUG("b4 add constraint");
        m_solver->addSafeSetConstraint(beta_coeffs, gamma_coeffs);
        DEBUG("af add constraint");
    }
    this->m_initialized = true;
}

template <std::size_t DIM>
BRY::Synthesizer<DIM>::Result BRY::PolyDynamicsSynthesizer<DIM>::synthesize(uint32_t time_horizon) {
    if (!this->m_initialized) {
        ERROR("Synthesizer not initialized (call initialize() before synthesize())");
        return typename Synthesizer<DIM>::Result{};
    }
    LPSolver::Result solver_result = m_solver->solve(time_horizon);

    typename Synthesizer<DIM>::Result result;
    result.p_safe = solver_result.p_safe;
    result.eta = solver_result.eta;
    result.gamma = solver_result.gamma;

    result.certificate.reset(new Polynomial<DIM, Basis::Bernstein>(solver_result.beta_values));
    return result;
}

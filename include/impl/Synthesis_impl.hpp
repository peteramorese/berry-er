#pragma once

#include "Synthesis.h"

#include "berry/BernsteinTransform.h"

/* Synthesizer */

template <std::size_t DIM>
void BRY::Synthesizer<DIM>::insertInitialSets(std::vector<HyperRectangle<DIM>>&& sets) {
    m_init_sets = std::move(sets);
}

template <std::size_t DIM>
void BRY::Synthesizer<DIM>::insertSafeSets(std::vector<HyperRectangle<DIM>>&& sets) {
    m_safe_sets = std::move(sets);
}

template <std::size_t DIM>
void BRY::Synthesizer<DIM>::insertUnsafeSets(std::vector<HyperRectangle<DIM>>&& sets) {
    m_unsafe_sets = std::move(sets);
}

/* PolyDynamicsSynthesizer */

template <std::size_t DIM>
BRY::PolyDynamicsSynthesizer<DIM>::PolyDynamicsSynthesizer(const std::shared_ptr<PolynomialDynamics<DIM>>& dynamics, const std::shared_ptr<Additive2ndMomentNoise<DIM>>& noise, bry_deg_t barrier_deg)
    : m_dynamics(dynamics)
    , m_solver(new LPSolver("SCIP", pow(barrier_deg + 1, DIM)))
    , m_barrier_deg(barrier_deg)
    , m_noise(noise)
{}

template <std::size_t DIM>
void BRY::PolyDynamicsSynthesizer<DIM>::initialize() {
    Eigen::MatrixXd Phi_m = BernsteinBasisTransform<DIM>::bernToPwrMatrix(m_barrier_deg);
    Eigen::MatrixXd Phi_inv_m = BernsteinBasisTransform<DIM>::pwrToBernMatrix(m_barrier_deg);
    DEBUG("Phi m: \n" << Phi_m);
    DEBUG("Phi inv m: \n" << Phi_inv_m);

    // Initial sets
    if (this->m_init_sets.empty())
        WARN("No initial sets were provided");
    //Eigen::VectorXd eta_coeffs = Phi_inv_m * Eigen::VectorXd::Ones(Phi_inv_m.cols());
    Eigen::VectorXd eta_coeffs = Eigen::VectorXd::Ones(Phi_inv_m.cols());
    for (const HyperRectangle<DIM>& set : this->m_init_sets) {
        DEBUG("Initial set transformation matrix \n" << set.transformationMatrix(m_barrier_deg));
        Eigen::MatrixXd beta_coeffs = -Phi_inv_m * set.transformationMatrix(m_barrier_deg) * Phi_m;
        DEBUG("Initial set beta coeffs: \n" << beta_coeffs);
        //DEBUG("Initial set eta coeffs: \n" << eta_coeffs);
        m_solver->addInitialSetConstraint(beta_coeffs, eta_coeffs);
    }

    // Unsafe sets
    if (this->m_unsafe_sets.empty())
        WARN("No unsafe sets were provided");
    Eigen::VectorXd lower_bound = Eigen::VectorXd::Ones(Phi_inv_m.cols());
    //Eigen::VectorXd lower_bound = Phi_inv_m * Eigen::VectorXd::Ones(Phi_inv_m.cols());
    //DEBUG("ones " << Phi_inv_m);
    //DEBUG("lower obund: " << lower_bound);
    for (const HyperRectangle<DIM>& set : this->m_unsafe_sets) {
        DEBUG("Set lower bound: " << set.lower_bounds.transpose());
        DEBUG("Unsafe set transformation matrix \n" << set.transformationMatrix(m_barrier_deg));
        Eigen::MatrixXd beta_coeffs = Phi_inv_m * set.transformationMatrix(m_barrier_deg) * Phi_m;
        DEBUG("Unsafe set beta coeffs: \n" << beta_coeffs);
        //DEBUG("Unsafe set eta coeffs: \n" << eta_coeffs);
        m_solver->addUnsafeSetConstraint(beta_coeffs, eta_coeffs);
    }

    // Safe sets
    if (this->m_safe_sets.empty())
        WARN("No safe sets were provided");
    Eigen::MatrixXd F_expec_Gamma = m_dynamics->dynamicsPowerMatrix(m_barrier_deg) * m_noise->additiveNoiseMatrix(m_barrier_deg);
    DEBUG("F \n" << m_dynamics->dynamicsPowerMatrix(m_barrier_deg));
    DEBUG("F expec Gamma: \n" << F_expec_Gamma);
    Eigen::MatrixXd Phi_inv_p = BernsteinBasisTransform<DIM>::pwrToBernMatrix(m_dynamics->composedDegree(m_barrier_deg));
    //Eigen::VectorXd gamma_coeffs = Phi_inv_p * Eigen::VectorXd::Ones(Phi_inv_p.cols());
    Eigen::VectorXd gamma_coeffs = Eigen::VectorXd::Ones(Phi_inv_p.cols());
    ASSERT(F_expec_Gamma.rows() == Phi_inv_p.cols(), "Dimension mismatch between F and Phi inv (p)");

    for (const HyperRectangle<DIM>& set : this->m_safe_sets) {
        Eigen::MatrixXd beta_coeffs = 
            -Phi_inv_p * (F_expec_Gamma - Eigen::MatrixXd::Identity(F_expec_Gamma.rows(), F_expec_Gamma.cols())) 
            * set.transformationMatrix(m_barrier_deg) * Phi_m;
        //DEBUG("Safe set beta coeffs: \n" << beta_coeffs);
        DEBUG("transformation matrix: \n" << set.transformationMatrix(m_barrier_deg));
        DEBUG("matrix: \n" << (F_expec_Gamma - Eigen::MatrixXd::Identity(F_expec_Gamma.rows(), F_expec_Gamma.cols())) * set.transformationMatrix(m_barrier_deg));
        m_solver->addSafeSetConstraint(beta_coeffs, gamma_coeffs);
    }
    this->m_initialized = true;
}

template <std::size_t DIM>
BRY::Synthesizer<DIM>::Result BRY::PolyDynamicsSynthesizer<DIM>::synthesize(uint32_t time_horizon) {
    if (!this->m_initialized) {
        ERROR("Synthesizer not initialized");
        return typename Synthesizer<DIM>::Result{};
    }
    LPSolver::Result solver_result = m_solver->solve(time_horizon);

    DEBUG("eta: " << solver_result.eta);
    DEBUG("gamma: " << solver_result.gamma);
    DEBUG("beta values: \n" << solver_result.beta_values);

    typename Synthesizer<DIM>::Result result;
    result.p_safe = solver_result.p_safe;
    result.certificate.reset(new Polynomial<DIM, Basis::Bernstein>(solver_result.beta_values));
    return result;
}

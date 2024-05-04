#pragma once

#include "Synthesis.h"

#include "berry/BernsteinTransform.h"

/* Synthesizer */

template <std::size_t DIM>
BRY::Synthesizer::insertInitialSets(std::vector<HyperRectangle<DIM>>&& sets) {
    m_init_sets = std::move(sets);
}

template <std::size_t DIM>
BRY::Synthesizer::insertSafeSets(std::vector<HyperRectangle<DIM>>&& sets) {
    m_safe_sets = std::move(sets);
}

template <std::size_t DIM>
BRY::Synthesizer::insertUnsafeSets(std::vector<HyperRectangle<DIM>>&& sets) {
    m_unsafe_sets = std::move(sets);
}

/* PolyDynamicsSynthesizer */

template <std::size_t DIM>
BRY::PolyDynamicsSynthesizer::PolyDynamicsSynthesizer(const std::shared_ptr<PolynomialDynamics<DIM>>& dynamics, const std::shared_ptr<Additive2ndMomentNoise<DIM>>& noise, bry_deg_t barrier_deg)
    : m_dynamics(dynamics)
    , m_solver(new LPSolver("SCIP", pow(barrier_deg + 1, DIM)))
    , m_barrier_deg(barrier_deg)
    , m_noise(noise)
{}

template <std::size_t DIM>
BRY::PolyDynamicsSynthesizer::initialize() {
    Eigen::MatrixXd Phi_m = BernsteinBasisTransform::getTfMatrix(m_barrier_deg);
    Eigen::MatrixXd Phi_inv_m = BernsteinBasisTransform::getInvTfMatrix(m_barrier_deg);

    // Initial sets
    Eigen::VectorXd eta_coeffs = Phi_inv_m * Eigen::VectorXd::Ones(Phi_inv_m.cols());
    for (const HyperRectangle<DIM>& set : m_init_sets) {
        Eigen::MatrixXd beta_coeffs = Phi_inv_m * set.getTransformationMatrix(m_barrier_deg) * Phi_m;
        m_solver->addInitialSetConstraint(-beta_coeffs, eta_coeffs);
    }

    // Unsafe sets
    Eigen::VectorXd lower_bound = Phi_inv_m * Eigen::VectorXd::Ones(Phi_inv_m.cols());
    for (const HyperRectangle<DIM>& set : m_unsafe_sets) {
        Eigen::MatrixXd beta_coeffs = Phi_inv_m * set.getTransformationMatrix(m_barrier_deg) * Phi_m;
        m_solver->addUnsafeSetConstraint(beta_coeffs, eta_coeffs);
    }

    // Safe sets
    Eigen::MatrixXd F_expec_Gamma = m_dynamics->dynamicsPowerMatrix(m_barrier_deg);
    Eigen::MatrixXd expec_Gamma = m_noise->additiveNoiseMatrix(m_barrier_deg);
    Eigen::MatrixXd Phi_inv_p = BernsteinBasisTransform::getInvTfMatrix(m_dynamics->composedDegree(m_barrier_deg));
    Eigen::VectorXd gamma_coeffs = Phi_inv_p * Eigen::VectorXd::Ones(Phi_inv_p.cols());
    ASSERT(F.rows() == Phi_inv_p.cols(), "Dimension mismatch between F and Phi inv (p)");
}

template <std::size_t DIM>
BRY::PolyDynamicsSynthesizer::synthesize(uint32_t time_horizon) {

}

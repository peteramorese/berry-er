#pragma once

#include "Synthesis.h"

#include "berry/BernsteinTransform.h"

/* SynthesisProblem */

template <std::size_t DIM>
void BRY::SynthesisProblem<DIM>::setWorkspace(const HyperRectangle<DIM>& workspace) {
    workspace_sets = {workspace};
}

template <std::size_t DIM>
std::unique_ptr<BRY::SynthesisProblem<DIM>> BRY::SynthesisProblem<DIM>::makeSubdividedProblem(uint32_t subdivision) {
    std::unique_ptr<BRY::SynthesisProblem<DIM>> subd_prob(new BRY::SynthesisProblem<DIM>());
    auto divide = [&] (std::vector<HyperRectangle<DIM>>& subd_sets, const std::vector<HyperRectangle<DIM>>& original_sets) {
        subd_sets.clear();
        subd_sets.reserve(original_sets.size() * pow(subdivision, DIM));
        for (const auto& set : original_sets) {
            std::vector<HyperRectangle<DIM>> subd_sets_to_insert = set.subdivide(subdivision);
            subd_sets.insert(subd_sets.end(), subd_sets_to_insert.begin(), subd_sets_to_insert.end());
        }
    };
    divide(subd_prob->workspace_sets, workspace_sets);
    divide(subd_prob->init_sets, init_sets);
    divide(subd_prob->safe_sets, safe_sets);
    divide(subd_prob->unsafe_sets, unsafe_sets);
    subd_prob->time_horizon = time_horizon;
    return subd_prob;
}

/* PolyDynamicsSynthesizer */

template <std::size_t DIM>
BRY::PolyDynamicsSynthesizer<DIM>::PolyDynamicsSynthesizer(const std::shared_ptr<PolynomialDynamics<DIM>>& dynamics, const std::shared_ptr<Additive2ndMomentNoise<DIM>>& noise, bry_int_t barrier_deg, const std::string& solver_id)
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
std::pair<BRY::Matrix, BRY::Vector> BRY::PolyDynamicsSynthesizer<DIM>::getConstraintMatrices(bry_int_t degree_increase, uint32_t subdivision) const {
    std::unique_ptr<SynthesisProblem<DIM>> subd_prob;
    SynthesisProblem<DIM>* problem;
    if (subdivision) {
        subd_prob = m_problem->makeSubdividedProblem(subdivision);
        problem = subd_prob.get();
    } else {
        problem = m_problem.get();
    }

    Matrix Phi_m = BernsteinBasisTransform<DIM>::pwrToBernMatrix(m_barrier_deg, degree_increase);

    bry_int_t n_cols = Phi_m.cols() + 2;

    // Workspace
    if (problem->workspace_sets.empty())
        WARN("No workspace was provided");
    Matrix ws_coeffs(problem->workspace_sets.size() * Phi_m.rows(), n_cols);
    Vector ws_lower_bound = Vector::Zero(ws_coeffs.rows());
    bry_int_t i = 0;
    for (const HyperRectangle<DIM>& set : problem->workspace_sets) {
        Matrix beta_coeffs = Phi_m * set.transformationMatrix(m_barrier_deg);
        Matrix coeffs(beta_coeffs.rows(), n_cols);
        coeffs << beta_coeffs, Vector::Zero(beta_coeffs.rows()), Vector::Zero(beta_coeffs.rows());
        ws_coeffs.block(Phi_m.rows() * i++, 0, Phi_m.rows(), n_cols) = coeffs;
    }

    // Initial sets
    if (problem->init_sets.empty())
        WARN("No initial sets were provided");
    Matrix init_coeffs(problem->init_sets.size() * Phi_m.rows(), n_cols);
    Vector init_lower_bound = Vector::Zero(init_coeffs.rows());
    i = 0;
    for (const HyperRectangle<DIM>& set : problem->init_sets) {
        Matrix beta_coeffs = -Phi_m * set.transformationMatrix(m_barrier_deg);
        Matrix coeffs(beta_coeffs.rows(), n_cols);
        coeffs << beta_coeffs, Vector::Ones(beta_coeffs.rows()), Vector::Zero(beta_coeffs.rows());
        init_coeffs.block(Phi_m.rows() * i++, 0, Phi_m.rows(), n_cols) = coeffs;
    }

    // Unsafe sets
    if (problem->unsafe_sets.empty())
        WARN("No unsafe sets were provided");
    Matrix unsafe_coeffs(problem->unsafe_sets.size() * Phi_m.rows(), n_cols);
    Vector unsafe_lower_bound = Vector::Ones(unsafe_coeffs.rows());
    i = 0;
    for (const HyperRectangle<DIM>& set : problem->unsafe_sets) {
        Matrix beta_coeffs = Phi_m * set.transformationMatrix(m_barrier_deg);
        Matrix coeffs(beta_coeffs.rows(), n_cols);
        coeffs << beta_coeffs, Vector::Zero(beta_coeffs.rows()), Vector::Zero(beta_coeffs.rows());
        unsafe_coeffs.block(Phi_m.rows() * i++, 0, Phi_m.rows(), n_cols) = coeffs;
    }

    // Safe sets
    if (problem->safe_sets.empty())
        WARN("No safe sets were provided");
    Matrix F_expec_Gamma = m_dynamics->dynamicsPowerMatrix(m_barrier_deg) * m_noise->additiveNoiseMatrix(m_barrier_deg);
    bry_int_t p = m_dynamics->composedDegree(m_barrier_deg);
    Matrix Phi_p = BernsteinBasisTransform<DIM>::pwrToBernMatrix(p, degree_increase);
    Vector gamma_coeffs = Vector::Ones(Phi_p.cols());
    ASSERT(F_expec_Gamma.rows() == Phi_p.cols(), "Dimension mismatch between F and Phi (p)");

    Matrix safe_coeffs(problem->safe_sets.size() * Phi_p.rows(), n_cols);
    Vector safe_lower_bound = Vector::Zero(safe_coeffs.rows());
    i = 0;
    for (const HyperRectangle<DIM>& set : problem->safe_sets) {
        Matrix beta_coeffs = 
            -Phi_p * set.transformationMatrix(p) * (F_expec_Gamma - Matrix::Identity(F_expec_Gamma.rows(), F_expec_Gamma.cols()));
        
        Matrix coeffs(beta_coeffs.rows(), n_cols);
        coeffs << beta_coeffs, Vector::Zero(beta_coeffs.rows()), Vector::Ones(beta_coeffs.rows());
        safe_coeffs.block(Phi_p.rows() * i++, 0, Phi_p.rows(), n_cols) = coeffs;
    }

    Matrix A(ws_coeffs.rows() + init_coeffs.rows() + unsafe_coeffs.rows() + safe_coeffs.rows(), n_cols);
    A << ws_coeffs, init_coeffs, unsafe_coeffs, safe_coeffs;

    Vector b(ws_lower_bound.size() + init_lower_bound.size() + unsafe_lower_bound.size() + safe_lower_bound.size());
    b << ws_lower_bound, init_lower_bound, unsafe_lower_bound, safe_lower_bound;
    return std::make_pair(std::move(A), std::move(b));
}

template <std::size_t DIM>
void BRY::PolyDynamicsSynthesizer<DIM>::initialize(bry_int_t degree_increase, uint32_t subdivision) {
    std::unique_ptr<SynthesisProblem<DIM>> subd_prob;
    SynthesisProblem<DIM>* problem;
    if (subdivision) {
        subd_prob = m_problem->makeSubdividedProblem(subdivision);
        problem = subd_prob.get();
    } else {
        problem = m_problem.get();
    }

    Matrix Phi_m = BernsteinBasisTransform<DIM>::pwrToBernMatrix(m_barrier_deg, degree_increase);
    Matrix Phi_inv_m = BernsteinBasisTransform<DIM>::bernToPwrMatrix(m_barrier_deg);

    // Workspace
    for (HyperRectangle<DIM> set : problem->workspace_sets) {
        Matrix beta_coeffs = Phi_m * set.transformationMatrix(m_barrier_deg) * Phi_inv_m;
        m_solver->addWorkspaceConstraint(beta_coeffs);
    }

    // Initial sets
    if (problem->init_sets.empty())
        WARN("No initial sets were provided");
    Vector eta_coeffs = Vector::Ones(Phi_m.rows());
    for (const HyperRectangle<DIM>& set : problem->init_sets) {
        Matrix beta_coeffs = -Phi_m * set.transformationMatrix(m_barrier_deg) * Phi_inv_m;
        m_solver->addInitialSetConstraint(beta_coeffs, eta_coeffs);
    }

    // Unsafe sets
    if (problem->unsafe_sets.empty())
        WARN("No unsafe sets were provided");
    Vector lower_bound = Vector::Ones(Phi_m.rows());
    for (const HyperRectangle<DIM>& set : problem->unsafe_sets) {
        Matrix beta_coeffs = Phi_m * set.transformationMatrix(m_barrier_deg) * Phi_inv_m;
        m_solver->addUnsafeSetConstraint(beta_coeffs, lower_bound);
    }

    // Safe sets
    if (problem->safe_sets.empty())
        WARN("No safe sets were provided");
    Matrix F_expec_Gamma = m_dynamics->dynamicsPowerMatrix(m_barrier_deg) * m_noise->additiveNoiseMatrix(m_barrier_deg);
    bry_int_t p = m_dynamics->composedDegree(m_barrier_deg);
    Matrix Phi_p = BernsteinBasisTransform<DIM>::pwrToBernMatrix(p, degree_increase);
    Vector gamma_coeffs = Vector::Ones(Phi_p.rows());
    ASSERT(F_expec_Gamma.rows() == Phi_p.cols(), "Dimension mismatch between F and Phi (p)");

    for (const HyperRectangle<DIM>& set : problem->safe_sets) {
        Matrix beta_coeffs = 
            -Phi_p * set.transformationMatrix(p) * (F_expec_Gamma - Matrix::Identity(F_expec_Gamma.rows(), F_expec_Gamma.cols())) * Phi_inv_m;

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

    solution.certificate.reset(new Polynomial<DIM, Basis::Bernstein>(solver_result.b_values));
    return solution;
}

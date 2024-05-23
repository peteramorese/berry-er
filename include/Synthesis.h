#pragma once

#include "HyperRectangle.h"
#include "Dynamics.h"
#include "Noise.h"
#include "LPSolver.h"

#include <vector>
#include <memory>

namespace BRY {

template <std::size_t DIM>
struct SynthesisProblem {
    HyperRectangle<DIM> workspace = HyperRectangle<DIM>();
    std::vector<HyperRectangle<DIM>> init_sets;
    std::vector<HyperRectangle<DIM>> safe_sets;
    std::vector<HyperRectangle<DIM>> unsafe_sets;
    uint32_t time_horizon;
};

template <std::size_t DIM>
struct SynthesisSolution {
    bool success = false;

    /// @brief Probability of safety
    bry_float_t p_safe = 0.0;

    /// @brief eta (init set constraint) and gamma (expected increase constraint)
    bry_float_t eta = 0.0; 
    bry_float_t gamma = 0.0;

    /// @brief Bernstein polynomial barrier certificate
    std::unique_ptr<Polynomial<DIM, Basis::Bernstein>> certificate = nullptr;
};

template <std::size_t DIM>
class PolyDynamicsSynthesizer {
    public:
        PolyDynamicsSynthesizer(const std::shared_ptr<PolynomialDynamics<DIM>>& dynamics, 
                                const std::shared_ptr<Additive2ndMomentNoise<DIM>>& noise, 
                                bry_int_t barrier_deg,
                                const std::string& solver_id = "SCIP");

        void setProblem(const std::shared_ptr<SynthesisProblem<DIM>>& problem);

        /// @brief Get the linear program matrices
        /// @return Pair: {A, b} where Ax >= b
        std::pair<Matrix, Vector> getConstraintMatrices(bry_int_t degree_increase = 0) const;

        /// @brief Set a time limit for the solver
        /// @param time_limit_ms Time limit in milliseconds
        void setTimeLimit(int64_t time_limit_ms);

        /// @brief Initialize the synthesizer
        /// @param degree_increase Decrease the conservativeness by increasing the number constraints
        void initialize(bry_int_t degree_increase = 0);

        SynthesisSolution<DIM> synthesize();

    private:
        std::shared_ptr<PolynomialDynamics<DIM>> m_dynamics;
        std::shared_ptr<Additive2ndMomentNoise<DIM>> m_noise;
        std::unique_ptr<LPSolver> m_solver;
        bry_int_t m_barrier_deg;

        bool m_initialized = false;

        std::shared_ptr<SynthesisProblem<DIM>> m_problem;
};

}

#include "impl/Synthesis_impl.hpp"
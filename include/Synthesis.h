#pragma once

#include "HyperRectangle.h"
#include "Dynamics.h"
#include "Noise.h"
#include "LPSolver.h"

#include <vector>
#include <memory>

namespace BRY {

template <std::size_t DIM>
class Synthesizer {
    public:
        struct Result {
            /// @brief Probability of safety
            bry_float_t p_safe;

            /// @brief Bernstein polynomial barrier certificate
            std::unique_ptr<Polynomial<DIM, Basis::Bernstein>> certificate;
        };

    public:

        void insertInitialSets(std::vector<HyperRectangle<DIM>>&& sets);
        void insertSafeSets(std::vector<HyperRectangle<DIM>>&& sets);
        void insertUnsafeSets(std::vector<HyperRectangle<DIM>>&& sets);
    
        virtual void initialize() = 0;
        virtual Result synthesize(uint32_t time_horizon) = 0;

    protected:
        Synthesizer() = default;

    protected:
        std::vector<HyperRectangle<DIM>> m_init_sets;
        std::vector<HyperRectangle<DIM>> m_safe_sets;
        std::vector<HyperRectangle<DIM>> m_unsafe_sets;
        bool m_initialized = false;
};

template <std::size_t DIM>
class PolyDynamicsSynthesizer : public Synthesizer<DIM> {
    public:
        PolyDynamicsSynthesizer(const std::shared_ptr<PolynomialDynamics<DIM>>& dynamics, 
                                const std::shared_ptr<Additive2ndMomentNoise<DIM>>& noise, 
                                bry_deg_t barrier_deg);

        virtual void initialize() override;
        virtual Synthesizer<DIM>::Result synthesize(uint32_t time_horizon) override;

    private:
        std::shared_ptr<PolynomialDynamics<DIM>> m_dynamics;
        std::shared_ptr<Additive2ndMomentNoise<DIM>> m_noise;
        std::unique_ptr<LPSolver> m_solver;
        bry_deg_t m_barrier_deg;
};

}

#include "impl/Synthesis_impl.hpp"
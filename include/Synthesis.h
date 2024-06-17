#pragma once

#include "HyperRectangle.h"
#include "Dynamics.h"
#include "Noise.h"
#include "LPSolver.h"

#include <vector>
#include <memory>

namespace BRY {

template <std::size_t DIM>
struct ConstraintMatrices {
    ConstraintMatrices(bry_int_t n_constraints, bry_int_t n_vars, bry_int_t barrier_deg, bool diag_deg);
    ConstraintMatrices() = delete;

    /// @brief Constraint matrix `A` in `Ax >= b` where columns represent variables [b_0, ..., b_m, eta, gamma]
    Matrix A;

    /// @brief Lower bound vector `b` in `Ax >= b`
    Vector b;

    /// @brief Degree definition of barrier
    BRY_INL bool diagDeg() const {return m_diag_deg;}

    void toDiagonalDegree();

    private:
        const bry_int_t m_barrier_deg;
        bool m_diag_deg;
};

template <std::size_t DIM>
struct PolyDynamicsProblem {

    /// @brief Dynamics
    std::shared_ptr<PolynomialDynamics<DIM>> dynamics;

    /// @brief Noise
    std::shared_ptr<Additive2ndMomentNoise<DIM>> noise;

    /* Set definitions */
    std::vector<HyperRectangle<DIM>> workspace_sets = {HyperRectangle<DIM>()};
    std::vector<HyperRectangle<DIM>> init_sets;
    std::vector<HyperRectangle<DIM>> safe_sets;
    std::vector<HyperRectangle<DIM>> unsafe_sets;

    /// @brief Number of time steps to verify the system for
    uint32_t time_horizon;

    /// @brief Degree of the barrier certificate
    bry_int_t barrier_deg = 3;

    /// @brief Increase the degree of the power-Bernstein conversion to reduce the conservativeness
    bry_int_t degree_increase = 0;

    /// @brief Use the diagonal degree definition of the certificate
    bool diag_deg;

    public:
        /// @brief Helper for setting the workspace from a single set
        /// @param workspace Hyper-rectangle workspace
        void setWorkspace(const HyperRectangle<DIM>& workspace);

        /// @brief Subdivide the problem for less conservativeness
        /// @param subdivision Integer number of subdivisions along one dimension
        void subdivide(uint32_t subdivision);

        /// @brief Compute the constraint matrices
        /// @return 
        const ConstraintMatrices getConstraintMatrices() const;
};

template <std::size_t DIM>
struct SynthesisSolution : public LPSolver::Result {
    SynthesisSolution(bool diag_deg) : m_diag_deg(diag_deg) {}

    /// @brief Degree definition of barrier
    BRY_INL bool diagDeg() const {return m_diag_deg;}

    void fromDiagonalDegree();

    private:
        bool m_diag_deg;
};

template <std::size_t DIM>
SynthesisSolution<DIM> synthesize(const PolyDynamicsProblem<DIM>& problem, bry_int_t barrier_deg, const std::string& solver_id = "clp");

}

#include "impl/Synthesis_impl.hpp"
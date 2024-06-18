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
    ConstraintMatrices(bry_int_t n_constraints, bry_int_t n_vars, bry_int_t barrier_deg_);
    ConstraintMatrices() = delete;

    /// @brief Constraint matrix `A` in `Ax >= b` where columns represent variables [b_0, ..., b_m, eta, gamma]
    Matrix A;

    /// @brief Lower bound vector `b` in `Ax >= b`
    Vector b;

    const bry_int_t barrier_deg;

    /// @brief Degree definition of barrier
    BRY_INL bool diagDeg() const {return m_diag_deg;}

    void toDiagonalDegree();

    private:
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
    uint32_t time_horizon = 10;

    /// @brief Degree of the barrier certificate
    bry_int_t barrier_deg = 3;

    /// @brief Increase the degree of the power-Bernstein conversion to reduce the conservativeness
    bry_int_t degree_increase = 0;

    /// @brief Use the diagonal degree definition of the certificate
    bool diag_deg = false;

    public:
        /// @brief Helper for setting the workspace from a single set
        /// @param workspace Hyper-rectangle workspace
        void setWorkspace(const HyperRectangle<DIM>& workspace);

        /// @brief Subdivide the problem for less conservativeness
        /// @param subdivision Integer number of subdivisions along one dimension
        void subdivide(uint32_t subdivision);

        /// @brief Compute the constraint matrices
        /// @return 
        const ConstraintMatrices<DIM> getConstraintMatrices() const;
};

template <std::size_t DIM>
struct SynthesisResult : public LPSolver::Result {
    SynthesisResult(bool diag_deg, bry_int_t barrier_deg_) : m_diag_deg(diag_deg), barrier_deg(barrier_deg_) {}

    /// @brief Degree definition of barrier
    BRY_INL bool diagDeg() const {return m_diag_deg;}

    void fromDiagonalDegree();

    const bry_int_t barrier_deg;

    private:
        bool m_diag_deg;
};

template <std::size_t DIM>
SynthesisResult<DIM> synthesize(const PolyDynamicsProblem<DIM>& problem, const std::string& solver_id = "clp");

template <std::size_t DIM>
SynthesisResult<DIM> synthesize(const ConstraintMatrices<DIM>& constraints, bry_int_t time_horizon, const std::string& solver_id = "clp");

void writeMatrixToFile(const Matrix& matrix, const std::string& filename);

}

#include "impl/Synthesis_impl.hpp"
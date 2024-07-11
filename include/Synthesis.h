#pragma once

#include "HyperRectangle.h"
#include "Dynamics.h"
#include "Noise.h"
#include "LPSolver.h"

#include <list>
#include <memory>

namespace BRY {

enum class ConstraintType {
    Workspace, 
    Init, 
    Safe, 
    Unsafe
};

struct ConstraintID {
    ConstraintType type;
    bry_int_t set_idx;
};

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

    /// @brief Convert the constraints to use diagonal degree
    void toDiagonalDegree();

    /// @brief Array with number elements equal to rows in `A` that identifies which set the constraint came from
    std::vector<ConstraintID> constraint_ids;

    private:
        bool m_diag_deg;
};

template <std::size_t DIM>
struct PolyDynamicsProblem {

    /// @brief Dynamics
    std::shared_ptr<PolynomialDynamics<DIM>> dynamics;

    /// @brief Noise
    std::shared_ptr<AdditiveGaussianNoise<DIM>> noise;

    /* Set definitions */
    std::list<HyperRectangle<DIM>> workspace_sets = {HyperRectangle<DIM>()};
    std::list<HyperRectangle<DIM>> init_sets;
    std::list<HyperRectangle<DIM>> safe_sets;
    std::list<HyperRectangle<DIM>> unsafe_sets;

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

        /// @brief Given a constraint ID, get the corresponding set that the constraint belongs to
        /// @param id ID of constraint
        /// @return Iterator to the set in the set list corresponding to the constraint type
        std::list<HyperRectangle<DIM>>::iterator lookupSetFromConstraint(const ConstraintID& id);
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

template <std::size_t DIM>
SynthesisResult<DIM> synthesizeAdaptive(const PolyDynamicsProblem<DIM>& problem, bry_int_t max_iter, bry_int_t subdiv_per_iter, const std::string& solver_id = "clp");

void writeMatrixToFile(const Matrix& matrix, const std::string& filename);

}

#include "impl/Synthesis_impl.hpp"
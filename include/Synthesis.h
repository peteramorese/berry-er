#pragma once

#include "HyperRectangle.h"
#include "Dynamics.h"
#include "Noise.h"
#include "LPSolver.h"
#include "MonomialFilter.h"

#include <list>
#include <memory>

namespace BRY {

enum ConstraintType {
    Workspace = 0, 
    Init = 1, 
    Unsafe = 2,
    Safe = 3, 
};

struct ConstraintID {
    ConstraintType type;
    bry_int_t set_idx;

    // For use as a map key
    BRY_INL bool operator<(const ConstraintID& other) const;
};

template <std::size_t DIM>
struct ConstraintMatrices {
    /// @brief Construct constraint matrices
    /// @param n_constraints Number of constriants
    /// @param n_vars Number of variables (number of monomial coefficients)
    /// @param barrier_deg_ Degree of the barrier
    ConstraintMatrices(bry_int_t n_constraints, bry_int_t n_vars, bry_int_t barrier_deg_);

    /// @brief Construct constraint matrices
    /// @param n_constraints Number of constriants
    /// @param n_vars Number of variables (number of monomial coefficients)
    /// @param barrier_deg_ Degree of the barrier
    /// @param filter Filter to apply to reduce the number of variables
    ConstraintMatrices(bry_int_t n_constraints, bry_int_t n_vars, bry_int_t barrier_deg_, const std::shared_ptr<MonomialFilter<DIM>>& filter_);

    ConstraintMatrices() = delete;

    /// @brief Constraint matrix `A` in `Ax >= b` where columns represent variables [b_0, ..., b_m, eta, gamma]
    Matrix A;

    /// @brief Lower bound vector `b` in `Ax >= b`
    Vector b;

    /// @brief Degree of barrier
    const bry_int_t barrier_deg;

    /// @brief Filter to apply to remove variables
    std::shared_ptr<MonomialFilter<DIM>> filter;

    /// @brief Array with number elements equal to rows in `A` that identifies which set the constraint came from
    std::vector<ConstraintID> constraint_ids;

    /// @brief Map from each constrained set to the transformation matrix. Each transformation does NOT include 
    /// the Bernstein basis conversion, i.e., power basis barrier to power basis polynomial constraint to lower-bound.
    /// Pointer will not own any object if the matrices were not stored upon construction (see 
    /// `PolyDynamicsProblem::getConstraintMatrices()`)
    std::unique_ptr<std::map<ConstraintID, Matrix>> transformation_matrices;

    /// @brief Determine if a filter has been applied to the matrices
    /// @return `true` if filter has been applied
    BRY_INL bool isFilterApplied() const {return m_filter_applied;}

    /// @brief Apply the current filter if one was supplied (otherwise do nothing)
    void applyFilter();

    /// @brief Compute the constraint robustness vector equal to `Av - b` where `v` is the solution vector. If
    /// the solution vector adheres to the constraints, the elements of the returned vector will be non-negative
    /// @param soln_vec Solution vector that optimizes the LP problem
    /// @return Robustness vector of size equal to number of rows in `A` where each element corresponds 
    /// to the robustness of a given constraint. 
    Vector computeRobustnessVec(const Vector& soln_vec) const;

    private:
        bool m_filter_applied;
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
    std::list<HyperRectangle<DIM>> unsafe_sets;
    std::list<HyperRectangle<DIM>> safe_sets;

    /// @brief Number of time steps to verify the system for
    uint32_t time_horizon = 10;

    /// @brief Degree of the barrier certificate
    bry_int_t barrier_deg = 3;

    /// @brief Increase the degree of the power-Bernstein conversion to reduce the conservativeness
    bry_int_t degree_increase = 0;

    /// @brief Monomial filter to apply to the constraints. If nullptr, no filter will be applied
    std::shared_ptr<MonomialFilter<DIM>> filter = nullptr;

    public:
        /// @brief Helper for setting the workspace from a single set
        /// @param workspace Hyper-rectangle workspace
        void setWorkspace(const HyperRectangle<DIM>& workspace);

        /// @brief Subdivide the problem for less conservativeness
        /// @param subdivision Integer number of subdivisions along one dimension
        void subdivide(uint32_t subdivision);

        /// @brief Get the total number of sets across workspace, init, unsafe, and safe.
        /// @return 
        BRY_INL bry_int_t numSets() const;

        /// @brief Compute the constraint matrices
        /// @param store_tf_matrices Returned object will cache the polynomial transformation matrices if true
        /// @return Constraint matrices object
        const ConstraintMatrices<DIM> getConstraintMatrices(bool store_tf_matrices = false) const;

        /// @brief Given a constraint ID, get the corresponding set that the constraint belongs to
        /// @param id ID of constraint
        /// @return Iterator to the set in the set list corresponding to the constraint type
        std::list<HyperRectangle<DIM>>::iterator lookupSetFromConstraint(const ConstraintID& id);
};

template <std::size_t DIM>
struct SynthesisResult : public LPSolver::Result {
    SynthesisResult() = default;
    SynthesisResult(bry_int_t barrier_deg_) : barrier_deg(barrier_deg_) {}
    SynthesisResult(bry_int_t barrier_deg_, const std::shared_ptr<MonomialFilter<DIM>>& filter_) : barrier_deg(barrier_deg_), filter(filter_) {}

    BRY_INL bool isFilterApplied() const {return (bool)filter;}

    /// @brief Remove the current filter if one was supplied (otherwise do nothing)
    void removeFilter();

    bry_int_t barrier_deg = 0;

    std::shared_ptr<MonomialFilter<DIM>> filter = nullptr;
};

template <std::size_t DIM>
SynthesisResult<DIM> synthesize(const PolyDynamicsProblem<DIM>& problem, const std::string& solver_id = "clp");

template <std::size_t DIM>
SynthesisResult<DIM> synthesize(const ConstraintMatrices<DIM>& constraints, bry_int_t time_horizon, const std::string& solver_id = "clp");

template <std::size_t DIM>
SynthesisResult<DIM> synthesizeAdaptive(PolyDynamicsProblem<DIM> problem, bry_int_t max_iter, bry_int_t subdiv_per_iter, const std::string& solver_id = "clp");

void writeMatrixToFile(const Matrix& matrix, const std::string& filename);

}

#include "impl/Synthesis_impl.hpp"
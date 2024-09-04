#pragma once

#include "HyperRectangle.h"
#include "HyperRectangle.h"
#include "Dynamics.h"
#include "Noise.h"

#include <map>

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

/// @brief Struct containing all of the set bounding information for a synthesis problem
template <std::size_t DIM>
struct SetDefinitions {
    
    std::multimap<ConstraintType, HyperRectangle<DIM>> sets;

    public:
        using Iterator = std::map<ConstraintType, HyperRectangle<DIM>>::iterator;
        using ConstIterator = std::map<ConstraintType, HyperRectangle<DIM>>::const_iterator;

    public:
        /// @brief Lex comparison forwarding for uniqueness checking and key ordering
        BRY_INL bool operator<(const SetDefinitions& other) const;

        std::pair<ConstIterator, ConstIterator> getSets(ConstraintType set_type) const;

        template <typename IT>
        void insertSets(ConstraintType set_type, IT begin, IT end);

        /// @brief Helper for setting the workspace from a single set
        /// @param workspace Hyper-rectangle workspace
        void setWorkspace(const HyperRectangle<DIM>& workspace);

        /// @brief Subdivide all sets 
        /// @param subdivision Integer number of subdivisions along one dimension
        void subdivide(uint32_t subdivision);
};

template <std::size_t DIM>
struct ConstraintMatrices {
    /// @brief Construct constraint matrices
    /// @param n_constraints Number of constriants
    /// @param n_vars Number of variables 
    /// @param barrier_deg_ Degree of the barrier
    ConstraintMatrices(bry_int_t n_constraints, bry_int_t n_vars, bry_int_t barrier_deg_);

    /// @brief Construct constraint matrices
    /// @param n_constraints Number of constriants
    /// @param n_vars Number of variables 
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

    /// @brief Filter applied to constraint matrices
    const std::shared_ptr<const MonomialFilter<DIM>> filter;

    /// @brief Array with number elements equal to rows in `A` that stores pointers to the set that created the constraint
    std::vector<typename SetDefinitions<DIM>::ConstIterator> constraint_sets;

    /// @brief Compute the constraint robustness vector equal to `Av - b` where `v` is the solution vector. If
    /// the solution vector adheres to the constraints, the elements of the returned vector will be non-negative
    /// @param soln_vec Solution vector that optimizes the LP problem
    /// @return Robustness vector of size equal to number of rows in `A` where each element corresponds 
    /// to the robustness of a given constraint. 
    Vector computeRobustnessVec(const Vector& soln_vec) const;
};

template <std::size_t DIM>
struct PolyDynamicsProblem : public SetDefinitions<DIM> {

    /// @brief Dynamics
    std::shared_ptr<PolynomialDynamics<DIM>> dynamics;

    /// @brief Noise
    std::shared_ptr<AdditiveGaussianNoise<DIM>> noise;

    /// @brief Number of time steps to verify the system for
    uint32_t time_horizon = 10;

    /// @brief Degree of the barrier certificate
    bry_int_t barrier_deg = 3;

    /// @brief Increase the degree of the power-Bernstein conversion to reduce the conservativeness
    bry_int_t degree_increase = 0;

    /// @brief Monomial filter to apply to the constraints. If nullptr, no filter will be applied
    std::shared_ptr<MonomialFilter<DIM>> filter = nullptr;

    public:
        /// @brief Compute the constraint matrices
        /// @return Constraint matrices object
        virtual const ConstraintMatrices<DIM> getConstraintMatrices();
    
        /// @brief Tighten the bounds on eta and gamma for a given barrier to get a better probability of safety
        /// @param result Result to edit (in place)
        /// @param eta_iterations Number of iterations to refine eta
        /// @param gamma_iterations Number of iterations to refine gamma
        void refineResult(LPSolver::Result& result, bry_int_t eta_iterations, bry_int_t gamma_iterations) const;

    protected:
        /// @brief Set `p`, `n_cols`, and `F_expec_Gamma_minus_I` before creating constraint matrices
        void initMatrixDefinitions();

        const Matrix& getPhim(bry_int_t bernstein_deg_incr) const;
        const Matrix& getPhip(bry_int_t bernstein_deg_incr) const;

        std::pair<Matrix, Vector> calculateSetConstraints(ConstraintType constraint_type, const HyperRectangle<DIM>& set) const;

    protected:
        bry_int_t m_p;
        bry_int_t m_n_cols;
        mutable std::map<bry_int_t, Matrix> m_Phi_m;
        mutable std::map<bry_int_t, Matrix> m_Phi_p;
        Matrix m_F_expec_Gamma_minus_I;

};

}

#include "impl/Problem_impl.hpp"
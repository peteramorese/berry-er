#pragma once

#include "Problem.h"

#include "berry/BernsteinTransform.h"

#include <fstream>
#include <iomanip>


bool BRY::ConstraintID::operator<(const ConstraintID& other) const {
    if (type == other.type) {
        return set_idx < other.set_idx;
    } else {
        return type < other.type;
    }
}

template <std::size_t DIM>
BRY::ConstraintMatrices<DIM>::ConstraintMatrices(bry_int_t n_constraints, bry_int_t n_vars, bry_int_t barrier_deg_)
    : A(n_constraints, n_vars)
    , b(n_constraints)
    , barrier_deg(barrier_deg_)
    , m_filter_applied(false)
{
    ASSERT(A.cols() == pow(barrier_deg + 1, DIM) + 2, "Number of vars does not match barrier degree + 2");
}

template <std::size_t DIM>
BRY::ConstraintMatrices<DIM>::ConstraintMatrices(bry_int_t n_constraints, bry_int_t n_vars, bry_int_t barrier_deg_, const std::shared_ptr<MonomialFilter<DIM>>& filter_)
    : A(n_constraints, n_vars)
    , b(n_constraints)
    , barrier_deg(barrier_deg_)
    , filter(filter_)
    , m_filter_applied(false)
{
    ASSERT(A.cols() == pow(barrier_deg + 1, DIM) + 2, "Number of vars does not match barrier degree + 2");
}

template <std::size_t DIM>
void BRY::ConstraintMatrices<DIM>::applyFilter() {
    if (m_filter_applied) {
        return;
    }

    bry_int_t n_coeffs = pow(barrier_deg + 1, DIM);

    bry_int_t removed_cols = 0;

    std::vector<bry_int_t> filter_flags(A.cols(), false);

    for (auto col_midx = mIdxW(DIM, barrier_deg + 1); !col_midx.last(); ++col_midx) {
        if (filter->remove(col_midx.begin())) {
            filter_flags[col_midx.inc().wrappedIdx()] = true;
            ++removed_cols;
        } 
    }

    Eigen::MatrixXd filtered_A(A.rows(), A.cols() - removed_cols);
    bry_int_t new_idx = 0;
    for (bry_int_t old_idx = 0; old_idx < A.cols(); ++old_idx) {
        if (!filter_flags[old_idx]) {
            filtered_A.col(new_idx++) = A.col(old_idx);
        }
    }

    A = filtered_A;

    m_filter_applied = true;
}

template <std::size_t DIM>
BRY::Vector BRY::ConstraintMatrices<DIM>::computeRobustnessVec(const Vector& soln_vec) const {
    return A * soln_vec - b;
}

template <std::size_t DIM>
std::pair<typename BRY::SetDefinitions<DIM>::ConstIterator, typename BRY::SetDefinitions<DIM>::ConstIterator> BRY::SetDefinitions<DIM>::getSets(ConstraintType set_type) const {
    return std::make_pair(sets.lower_bound(set_type), sets.upper_bound(set_type));
}

template <std::size_t DIM>
template <typename IT>
void BRY::SetDefinitions<DIM>::insertSets(ConstraintType set_type, IT begin, IT end) {
    for (IT it = begin; it != end; ++it) {
        sets.insert(std::make_pair(set_type, *it));
    }
}

template <std::size_t DIM>
void BRY::SetDefinitions<DIM>::setWorkspace(const HyperRectangle<DIM>& workspace) {
    sets.insert(std::make_pair(ConstraintType::Workspace, workspace));
}

template <std::size_t DIM>
void BRY::SetDefinitions<DIM>::subdivide(uint32_t subdivision) {
    if (subdivision < 2) {
        WARN("Subdivision is less than 2 (no effect)");
        return;
    }
    const std::multimap<ConstraintType, HyperRectangle<DIM>> original_sets = sets;
    sets.clear();
    for (const auto&[type, original_set] : original_sets) {
        std::vector<HyperRectangle<DIM>> subd_sets_to_insert = original_set.subdivide(subdivision);
        insertSets(type, subd_sets_to_insert.begin(), subd_sets_to_insert.end());
    }
}

template <std::size_t DIM>
const BRY::ConstraintMatrices<DIM> BRY::PolyDynamicsProblem<DIM>::getConstraintMatrices(bool store_tf_matrices) const {
    INFO("Creating constraint matrices");

    // Degree of the composed polynomial
    bry_int_t p = dynamics->composedDegree(barrier_deg);
    // Product F * E[Gamma]
    Matrix F_expec_Gamma = dynamics->dynamicsPowerMatrix(barrier_deg) * noise->additiveNoiseMatrix(barrier_deg);
    // Degree lift transform for the subtraction of F * E[Gamma] - B
    Matrix deg_lift_tf = makeDegreeChangeTransform<DIM>(barrier_deg, p);

    // Map containing all of the Phi_m transformation matrices for each bernstein conversion degree increase to prevent duplicate Phi matrices
    std::map<bry_int_t, Matrix> Phi_m;
    std::map<bry_int_t, Matrix> Phi_p;

    // Keep track of the number of constraints for each type
    std::array<bry_int_t, 4> n_constraints = makeUniformArray<bry_int_t, 4>(0);
    //bry_int_t n_ws_constraints = 0, n_init_constraints = 0, n_unsafe_constraints = 0, n_safe_constraints = 0;

    for (const auto&[set_type, set] : this->sets) {
        // Safe set constraints use the Phi_p container, all other constraints use Phi_m
        std::map<bry_int_t, Matrix>& Phi_matrix_container = set_type != ConstraintType::Safe ? Phi_m : Phi_p;

        bry_int_t Phi_deg = set_type != ConstraintType::Safe ? barrier_deg : p;

        /// Check if there is a Phi with the corresponding degree increase 
        auto it = Phi_matrix_container.find(set.bernstein_deg_incr);
        if (it == Phi_matrix_container.end()) {
            it = Phi_matrix_container.emplace(set.bernstein_deg_incr, BernsteinBasisTransform<DIM>::pwrToBernMatrix(Phi_deg, set.bernstein_deg_incr)).first;
        } 
        // Count the number of constraints for matrix allocation later
        n_constraints[set_type] += it->second.rows();
    }

    // Phi_m is guaranteed to atleast have one matrix in it, and all of the matrices must have the same number of columns, so look it up there
    bry_int_t n_cols = Phi_m.begin()->second.cols() + 2;

    BRY::ConstraintMatrices<DIM> constraint_matrices(std::accumulate(n_constraints.begin(), n_constraints.end(), 0), n_cols, barrier_deg);

    bry_int_t constraint_idx = 0;
    //set_idx = 0;
    for (const auto&[set_type, set] : this->sets) {
        if (set_type != ConstraintType::Safe) {
            // Eta coeffs are 1 if the set is an initial set, otherwise they are zero
            bry_float_t eta_coeff = static_cast<bry_float_t>(set_type == ConstraintType::Init);
            // Lower bound is zero unless the set is of type unsafe, in which case the lb is 1
            bry_float_t lb_coeff = static_cast<bry_float_t>(set_type == ConstraintType::Unsafe);
            // If the set is an initial set, then negate the b coefficients
            bry_float_t b_coeff_multiplier = (set_type == ConstraintType::Init) ? -1.0 : 1.0;

            Matrix tf = b_coeff_multiplier * set.transformationMatrix(barrier_deg);
            Matrix b_coeffs = Phi_m.at(set.bernstein_deg_incr) * tf;

            Matrix A_mat_vals(b_coeffs.rows(), n_cols);

            //            b         eta                                           gamma
            A_mat_vals << b_coeffs, Vector::Constant(b_coeffs.rows(), eta_coeff), Vector::Zero(b_coeffs.rows());
            constraint_matrices.A.block(constraint_idx, 0, A_mat_vals.rows(), A_mat_vals.cols()) = A_mat_vals;
            constraint_matrices.b.segment(constraint_idx, A_mat_vals.rows()) = Vector::Constant(A_mat_vals.rows(), lb_coeff);
            constraint_idx += A_mat_vals.rows();
        } else {
            Matrix tf = -set.transformationMatrix(p) * (F_expec_Gamma - deg_lift_tf);
            Matrix b_coeffs = Phi_p.at(set.bernstein_deg_incr) * tf;

            Matrix A_mat_vals(b_coeffs.rows(), n_cols);

            //            b         eta                            gamma
            A_mat_vals << b_coeffs, Vector::Zero(b_coeffs.rows()), Vector::Ones(b_coeffs.rows());
            constraint_matrices.A.block(constraint_idx, 0, A_mat_vals.rows(), A_mat_vals.cols()) = A_mat_vals;
            constraint_matrices.b.segment(constraint_idx, A_mat_vals.rows()) = Vector::Zero(A_mat_vals.rows());
            constraint_idx += A_mat_vals.rows();
        }
    }


    constraint_matrices.filter = filter;
    if (filter) {
        INFO("Applying filter...");
        constraint_matrices.applyFilter();
        INFO("Done!");
    }
    INFO("Created constraint matrices");
    return constraint_matrices;
}

template <std::size_t DIM>
std::list<BRY::HyperRectangle<DIM>>::iterator BRY::PolyDynamicsProblem<DIM>::lookupSetFromConstraint(const ConstraintID& id) {
    switch (id.type) {
        case ConstraintType::Workspace: {
            #ifdef BRY_ENABLE_BOUNDS_CHECK
                ASSERT(id.set_idx < this->workspace_sets.size(), "Set idx out of bounds (workspace sets)");
            #endif
            return std::next(this->workspace_sets.begin(), id.set_idx);
        }
        case ConstraintType::Init: {
            #ifdef BRY_ENABLE_BOUNDS_CHECK
                ASSERT(id.set_idx < this->init_sets.size(), "Set idx out of bounds (init sets)");
            #endif
            return std::next(this->init_sets.begin(), id.set_idx);
        }
        case ConstraintType::Unsafe: {
            #ifdef BRY_ENABLE_BOUNDS_CHECK
                ASSERT(id.set_idx < this->unsafe_sets.size(), "Set idx out of bounds (unsafe sets)");
            #endif
            return std::next(this->unsafe_sets.begin(), id.set_idx);
        }
        case ConstraintType::Safe: {
            #ifdef BRY_ENABLE_BOUNDS_CHECK
                ASSERT(id.set_idx < this->safe_sets.size(), "Set idx out of bounds (safe sets)");
            #endif
            return std::next(this->safe_sets.begin(), id.set_idx);
        }
    }
    throw std::invalid_argument("ID is invalid");
}
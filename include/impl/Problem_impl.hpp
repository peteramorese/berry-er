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
    , constraint_ids(n_constraints)
    , m_filter_applied(false)
{
    ASSERT(A.cols() == pow(barrier_deg + 1, DIM) + 2, "Number of vars does not match barrier degree + 2");
}

template <std::size_t DIM>
BRY::ConstraintMatrices<DIM>::ConstraintMatrices(bry_int_t n_constraints, bry_int_t n_vars, bry_int_t barrier_deg_, const std::shared_ptr<MonomialFilter<DIM>>& filter_)
    : A(n_constraints, n_vars)
    , b(n_constraints)
    , barrier_deg(barrier_deg_)
    , constraint_ids(n_constraints)
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
void BRY::PolyDynamicsProblem<DIM>::setWorkspace(const HyperRectangle<DIM>& workspace) {
    workspace_sets = {workspace};
}

template <std::size_t DIM>
void BRY::PolyDynamicsProblem<DIM>::subdivide(uint32_t subdivision) {
    if (subdivision < 2) {
        WARN("Subdivision is less than 2 (no effect)");
        return;
    }
    auto divide = [&] (std::list<HyperRectangle<DIM>>& subd_sets) {
        const std::list<HyperRectangle<DIM>> original_sets = subd_sets;
        subd_sets.clear();
        for (const auto& set : original_sets) {
            std::vector<HyperRectangle<DIM>> subd_sets_to_insert = set.subdivide(subdivision);
            subd_sets.insert(subd_sets.end(), subd_sets_to_insert.begin(), subd_sets_to_insert.end());
        }
    };
    divide(workspace_sets);
    divide(init_sets);
    divide(safe_sets);
    divide(unsafe_sets);
}

template <std::size_t DIM>
BRY::bry_int_t BRY::PolyDynamicsProblem<DIM>::numSets() const {
    return workspace_sets.size() + init_sets.size() + unsafe_sets.size() + safe_sets.size();
}

template <std::size_t DIM>
const BRY::ConstraintMatrices<DIM> BRY::PolyDynamicsProblem<DIM>::getConstraintMatrices(bool store_tf_matrices) const {
    INFO("Creating constraint matrices");

    Matrix Phi_m = BernsteinBasisTransform<DIM>::pwrToBernMatrix(barrier_deg, degree_increase);

    std::unique_ptr<std::map<ConstraintID, Matrix>> cached_constraint_matrices(store_tf_matrices ? new std::map<ConstraintID, Matrix>() : nullptr);

    bry_int_t n_cols = Phi_m.cols() + 2;

    // Workspace
    if (workspace_sets.empty())
        WARN("No workspace was provided");
    Matrix ws_coeffs(workspace_sets.size() * Phi_m.rows(), n_cols);
    Vector ws_lower_bound = Vector::Zero(ws_coeffs.rows());
    bry_int_t i = 0;
    for (const HyperRectangle<DIM>& set : workspace_sets) {
        Matrix tf = set.transformationMatrix(barrier_deg);
        Matrix beta_coeffs = Phi_m * tf;

        if (store_tf_matrices) {
            bool inserted = cached_constraint_matrices->emplace(ConstraintID{ConstraintType::Workspace, i}, std::move(tf)).second;
            ASSERT(inserted, "Duplicate constraint ID found");
        }

        Matrix coeffs(beta_coeffs.rows(), n_cols);
        coeffs << beta_coeffs, Vector::Zero(beta_coeffs.rows()), Vector::Zero(beta_coeffs.rows());
        ws_coeffs.block(Phi_m.rows() * i++, 0, Phi_m.rows(), n_cols) = coeffs;
    }
    INFO("Workspace constraints done. Computing initial set constraints...");

    // Initial sets
    if (init_sets.empty())
        WARN("No initial sets were provided");
    Matrix init_coeffs(init_sets.size() * Phi_m.rows(), n_cols);
    Vector init_lower_bound = Vector::Zero(init_coeffs.rows());
    i = 0;
    for (const HyperRectangle<DIM>& set : init_sets) {
        Matrix tf = -set.transformationMatrix(barrier_deg);
        Matrix beta_coeffs = Phi_m * tf;

        if (store_tf_matrices) {
            bool inserted = cached_constraint_matrices->emplace(ConstraintID{ConstraintType::Init, i}, std::move(tf)).second;
            ASSERT(inserted, "Duplicate constraint ID found");
        }

        Matrix coeffs(beta_coeffs.rows(), n_cols);
        coeffs << beta_coeffs, Vector::Ones(beta_coeffs.rows()), Vector::Zero(beta_coeffs.rows());
        init_coeffs.block(Phi_m.rows() * i++, 0, Phi_m.rows(), n_cols) = coeffs;
    }
    INFO("Initial sets done. Computing unsafe set constraints...");

    // Unsafe sets
    if (unsafe_sets.empty())
        WARN("No unsafe sets were provided");
    Matrix unsafe_coeffs(unsafe_sets.size() * Phi_m.rows(), n_cols);
    Vector unsafe_lower_bound = Vector::Ones(unsafe_coeffs.rows());
    i = 0;
    for (const HyperRectangle<DIM>& set : unsafe_sets) {
        Matrix tf = set.transformationMatrix(barrier_deg);
        Matrix beta_coeffs = Phi_m * tf;

        if (store_tf_matrices) {
            bool inserted = cached_constraint_matrices->emplace(ConstraintID{ConstraintType::Unsafe, i}, std::move(tf)).second;
            ASSERT(inserted, "Duplicate constraint ID found");
        }

        Matrix coeffs(beta_coeffs.rows(), n_cols);
        coeffs << beta_coeffs, Vector::Zero(beta_coeffs.rows()), Vector::Zero(beta_coeffs.rows());
        unsafe_coeffs.block(Phi_m.rows() * i++, 0, Phi_m.rows(), n_cols) = coeffs;
    }
    INFO("Unsafe sets done. Computing safe set constraints...");

    // Safe sets
    if (safe_sets.empty())
        WARN("No safe sets were provided");
    //DEBUG("Dynamics power matrix: \n" << dynamics->dynamicsPowerMatrix(barrier_deg));
    //DEBUG("Noise matrix: \n" << noise->additiveNoiseMatrix(barrier_deg));
    DEBUG("b4 f expec gamma compute");
    Matrix F_expec_Gamma = dynamics->dynamicsPowerMatrix(barrier_deg) * noise->additiveNoiseMatrix(barrier_deg);
    DEBUG("af f expec gamma compute");
    bry_int_t p = dynamics->composedDegree(barrier_deg);
    DEBUG("b4 phi p");
    Matrix Phi_p = BernsteinBasisTransform<DIM>::pwrToBernMatrix(p, degree_increase);
    DEBUG("af phi p (size: " << Phi_p.size());
    Vector gamma_coeffs = Vector::Ones(Phi_p.cols());
    ASSERT(F_expec_Gamma.rows() == Phi_p.cols(), "Dimension mismatch between F and Phi (p)");

    Matrix safe_coeffs(safe_sets.size() * Phi_p.rows(), n_cols);
    Vector safe_lower_bound = Vector::Zero(safe_coeffs.rows());
    i = 0;

    Matrix deg_lift_tf = makeDegreeChangeTransform<DIM>(barrier_deg, p);

    for (const HyperRectangle<DIM>& set : safe_sets) {
        DEBUG("b4 tf");
        Matrix tf = -set.transformationMatrix(p) * (F_expec_Gamma - deg_lift_tf);
        DEBUG("af tf");
        Matrix beta_coeffs = Phi_p * tf;

        if (store_tf_matrices) {
            bool inserted = cached_constraint_matrices->emplace(ConstraintID{ConstraintType::Safe, i}, std::move(tf)).second;
            ASSERT(inserted, "Duplicate constraint ID found");
        }
        
        Matrix coeffs(beta_coeffs.rows(), n_cols);
        coeffs << beta_coeffs, Vector::Zero(beta_coeffs.rows()), Vector::Ones(beta_coeffs.rows());
        safe_coeffs.block(Phi_p.rows() * i++, 0, Phi_p.rows(), n_cols) = coeffs;
    }
    INFO("Safe sets done.");

    BRY::ConstraintMatrices<DIM> constraint_matrices(ws_coeffs.rows() + init_coeffs.rows() + unsafe_coeffs.rows() + safe_coeffs.rows(), n_cols, barrier_deg);

    constraint_matrices.A << ws_coeffs, init_coeffs, unsafe_coeffs, safe_coeffs;
    constraint_matrices.b << ws_lower_bound, init_lower_bound, unsafe_lower_bound, safe_lower_bound;

    // Fill the constraint IDs
    auto it = constraint_matrices.constraint_ids.begin();
    for (bry_int_t set_i = 0; set_i < workspace_sets.size(); ++set_i) {
        std::fill(it, std::next(it, Phi_m.rows()), ConstraintID{ConstraintType::Workspace, set_i});
        std::advance(it, Phi_m.rows());
    }
    for (bry_int_t set_i = 0; set_i < init_sets.size(); ++set_i) {
        std::fill(it, std::next(it, Phi_m.rows()), ConstraintID{ConstraintType::Init, set_i});
        std::advance(it, Phi_m.rows());
    }
    for (bry_int_t set_i = 0; set_i < unsafe_sets.size(); ++set_i) {
        std::fill(it, std::next(it, Phi_m.rows()), ConstraintID{ConstraintType::Unsafe, set_i});
        std::advance(it, Phi_m.rows());
    }
    for (bry_int_t set_i = 0; set_i < safe_sets.size(); ++set_i) {
        std::fill(it, std::next(it, Phi_p.rows()), ConstraintID{ConstraintType::Safe, set_i}); // Advance by Phi_p rows instead of Phi_m
        std::advance(it, Phi_p.rows());
    }

    constraint_matrices.transformation_matrices = std::move(cached_constraint_matrices);
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
                ASSERT(id.set_idx < workspace_sets.size(), "Set idx out of bounds (workspace sets)");
            #endif
            return std::next(workspace_sets.begin(), id.set_idx);
        }
        case ConstraintType::Init: {
            #ifdef BRY_ENABLE_BOUNDS_CHECK
                ASSERT(id.set_idx < init_sets.size(), "Set idx out of bounds (init sets)");
            #endif
            return std::next(init_sets.begin(), id.set_idx);
        }
        case ConstraintType::Unsafe: {
            #ifdef BRY_ENABLE_BOUNDS_CHECK
                ASSERT(id.set_idx < unsafe_sets.size(), "Set idx out of bounds (unsafe sets)");
            #endif
            return std::next(unsafe_sets.begin(), id.set_idx);
        }
        case ConstraintType::Safe: {
            #ifdef BRY_ENABLE_BOUNDS_CHECK
                ASSERT(id.set_idx < safe_sets.size(), "Set idx out of bounds (safe sets)");
            #endif
            return std::next(safe_sets.begin(), id.set_idx);
        }
    }
    throw std::invalid_argument("ID is invalid");
}
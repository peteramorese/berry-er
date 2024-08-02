#pragma once

#include "AdaptiveProblem.h"

template <std::size_t DIM>
const BRY::ConstraintMatrices<DIM> BRY::AdaptiveProblem<DIM>::getConstraintMatrices(bool store_tf_matrices) const override {

    if (!existing_result) {
        return PolyDynamicsProblem<DIM>::
    }
    INFO("Creating constraint matrices");


    std::unique_ptr<std::map<ConstraintID, Matrix>> cached_constraint_matrices(store_tf_matrices ? new std::map<ConstraintID, Matrix>() : nullptr);


    // Things fixed before the algorithm
    bry_int_t n_cols = Phi_m.cols() + 2;
    Matrix F_expec_Gamma = dynamics->dynamicsPowerMatrix(barrier_deg) * noise->additiveNoiseMatrix(barrier_deg);
    bry_int_t p = dynamics->composedDegree(barrier_deg);

    Matrix Phi_m = BernsteinBasisTransform<DIM>::pwrToBernMatrix(barrier_deg, degree_increase);
    Matrix Phi_p = BernsteinBasisTransform<DIM>::pwrToBernMatrix(p, degree_increase);

    ASSERT(F_expec_Gamma.rows() == Phi_p.cols(), "Dimension mismatch between F and Phi (p)");

    Vector gamma_coeffs = Vector::Ones(Phi_p.cols());

    auto bry_int_t nConstraints = [&]() {
        return Phi_m.rows() * (workspace_sets.size() + init_sets.size() + unsafe_sets.size()) + Phi_p.rows() * (safe_sets.size());
    };

    while (nConstraints() < max_constraints) {}
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

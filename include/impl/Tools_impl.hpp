#pragma once

#include "Tools.h"

template <std::size_t DIM>
std::list<BRY::HyperRectangle<DIM>> BRY::makeRectBoundary(const HyperRectangle<DIM>& workspace, bry_float_t boundary_width, bry_int_t degree_increase) {
    std::list<BRY::HyperRectangle<DIM>> boundary_sets;
    for (bry_int_t d = 0; d < DIM; ++d) {
        Eigen::Vector<bry_float_t, DIM> boundary_width_vec;
        boundary_width_vec.setConstant(boundary_width);

        // Inflate the workspace by the boundary width
        HyperRectangle<DIM> h_lower = workspace;
        h_lower.bernstein_deg_incr = degree_increase;
        h_lower.lower_bounds -= boundary_width_vec;
        h_lower.upper_bounds += boundary_width_vec;

        // Copy over inflated workspace before editing the dim-specifc bounds
        HyperRectangle<DIM> h_upper = h_lower;

        // Chop off the part the overlaps the workspace
        h_lower.upper_bounds[d] = workspace.lower_bounds[d];
        h_upper.lower_bounds[d] = workspace.upper_bounds[d];

        // Add both sets to the boundary
        boundary_sets.push_back(std::move(h_lower));
        boundary_sets.push_back(std::move(h_upper));
    }
    return boundary_sets;
}

template <std::size_t DIM>
static void BRY::removeFilterFromPVector(Vector& polynomial_vector, bry_int_t barrier_deg, const MonomialFilter<DIM>* filter) {
    if (!filter) {
        return;
    }

    bry_int_t n_coeffs = pow(barrier_deg + 1, DIM);
    ASSERT(polynomial_vector.size() != n_coeffs, "Supplied coefficient vector is already in square-degree form (no filter is applied)");

    Vector sq_coeffs = Vector::Zero(n_coeffs);
    bry_int_t i = 0;

    for (auto col_midx = mIdxW(DIM, barrier_deg + 1); !col_midx.last(); ++col_midx) {
        if (!filter->remove(col_midx.begin())) {
            sq_coeffs(col_midx.inc().wrappedIdx()) = polynomial_vector(i++);
        }
    }

    polynomial_vector = sq_coeffs;
}
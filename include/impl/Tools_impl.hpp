#pragma once

#include "Tools.h"

template <std::size_t DIM>
std::list<BRY::HyperRectangle<DIM>> BRY::makeRectBoundary(const HyperRectangle<DIM>& workspace, bry_float_t boundary_width) {
    std::list<BRY::HyperRectangle<DIM>> boundary_sets;
    for (bry_int_t d = 0; d < DIM; ++d) {
        Eigen::Vector<bry_float_t, DIM> boundary_width_vec;
        boundary_width_vec.setConstant(boundary_width);

        // Inflate the workspace by the boundary width
        HyperRectangle<DIM> h_lower = workspace;
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
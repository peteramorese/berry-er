#pragma once

#include "HyperRectangle.h"

#include "berry/Options.h"

#include <Eigen/Core>

template <std::size_t DIM>
Eigen::Vector<BRY::bry_float_t, DIM> BRY::HyperRectangle<DIM>::scaleToUnit() const {
    Eigen::Vector<BRY::bry_float_t, DIM> diff = upper_bounds - lower_bounds;
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT((diff.array() > 0.0).all(), "Upper bounds are not greater than lower bounds");
#endif
    return diff.array().pow(-1.0);
}

template <std::size_t DIM>
Eigen::Vector<BRY::bry_float_t, DIM> BRY::HyperRectangle<DIM>::translationToUnit() const {
    return -lower_bounds.cwiseProduct(scaleToUnit());
}

template <std::size_t DIM>
Eigen::Vector<BRY::bry_float_t, DIM> BRY::HyperRectangle<DIM>::scaleFromUnit() const {
    Eigen::Vector<BRY::bry_float_t, DIM> diff = upper_bounds - lower_bounds;
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT((diff.array() > 0.0).all(), "Upper bounds are not greater than lower bounds");
#endif
    return diff;
}

template <std::size_t DIM>
Eigen::Vector<BRY::bry_float_t, DIM> BRY::HyperRectangle<DIM>::translationFromUnit() const {
    return lower_bounds;
}
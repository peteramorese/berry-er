#pragma once

#include "berry/Options.h"

namespace BRY {

template <std::size_t DIM>
struct HyperRectangle {
    public:
        Eigen::Vector<bry_float_t, DIM> lower_bounds;
        Eigen::Vector<bry_float_t, DIM> upper_bounds;

    public:
        BRY_INL Eigen::Vector<bry_float_t, DIM> scaleToUnit() const;
        BRY_INL Eigen::Vector<bry_float_t, DIM> translationToUnit() const;
        BRY_INL Eigen::Vector<bry_float_t, DIM> scaleFromUnit() const;
        BRY_INL Eigen::Vector<bry_float_t, DIM> translationFromUnit() const;
};

}

#include "impl/HyperRectangle_impl.hpp"
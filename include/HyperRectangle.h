#pragma once

#include "berry/Options.h"

namespace BRY {

template <std::size_t DIM>
struct HyperRectangle {
    public:
        Eigen::Vector<bry_float_t, DIM> lower_bounds = Eigen::Vector<bry_float_t, DIM>::Constant(0.0);
        Eigen::Vector<bry_float_t, DIM> upper_bounds = Eigen::Vector<bry_float_t, DIM>::Constant(1.0);

    public:
        Eigen::MatrixXd transformationMatrix(bry_deg_t m) const;

        BRY_INL Eigen::Vector<bry_float_t, DIM> scaleToUnit() const;
        BRY_INL Eigen::Vector<bry_float_t, DIM> translationToUnit() const;
        BRY_INL Eigen::Vector<bry_float_t, DIM> scaleFromUnit() const;
        BRY_INL Eigen::Vector<bry_float_t, DIM> translationFromUnit() const;
};

}

#include "impl/HyperRectangle_impl.hpp"
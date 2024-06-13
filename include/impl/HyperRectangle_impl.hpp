#pragma once

#include "HyperRectangle.h"

#include "berry/Options.h"
#include "berry/MultiIndex.h"

#include <Eigen/Core>

template <std::size_t DIM>
BRY::Matrix BRY::HyperRectangle<DIM>::transformationMatrix(bry_int_t m) const {
    BRY::bry_int_t m_monoms = pow(m + 1, DIM);

    Matrix T(m_monoms, m_monoms);
    T.setZero();

    Eigen::Vector<bry_float_t, DIM> scale = scaleFromUnit();
    Eigen::Vector<bry_float_t, DIM> translation = translationFromUnit();

    for (auto row_midx = mIdxW(DIM, m + 1); !row_midx.last(); ++row_midx) {
        std::vector<bry_int_t> index_bounds(row_midx.size());
        for (std::size_t d = 0; d < DIM; ++d) {
            index_bounds[d] = m + 1 - row_midx[d];
        }

        for (auto col_midx = mIdxBEW(index_bounds, m + 1); !col_midx.last(); ++col_midx) {
            bry_float_t element = 1.0;
            for (std::size_t j = 0; j < DIM; ++j) {
                element *= binom(row_midx[j] + col_midx[j], row_midx[j]);
                element *= std::pow(scale[j], row_midx[j]);
                element *= std::pow(translation[j], col_midx[j]);
            }

            T(row_midx.inc().wrappedIdx(), row_midx.inc().wrappedIdx() + col_midx.inc().wrappedIdx()) = element;
        }
    }
    return T;
}

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

template <std::size_t DIM>
std::vector<BRY::HyperRectangle<DIM>> BRY::HyperRectangle<DIM>::subdivide(uint32_t integer_division) const {
    std::vector<BRY::HyperRectangle<DIM>> subsets;
    subsets.reserve(pow(integer_division, DIM));
    Eigen::Vector<bry_float_t, DIM> diff = (upper_bounds - lower_bounds) / static_cast<float>(integer_division);
    for (auto midx = mIdx(DIM, integer_division); !midx.last(); ++midx) {
        Eigen::Vector<bry_float_t, DIM> shift = Eigen::Vector<bry_float_t, DIM>::Zero();
        for (std::size_t d = 0; d < DIM; ++d) {
            shift[d] = midx[d] * diff[d];
        }
        BRY::HyperRectangle<DIM> subset;
        subset.lower_bounds = lower_bounds + shift;
        subset.upper_bounds = subset.lower_bounds + diff;
        subsets.push_back(std::move(subset));
    }
    return subsets;
}
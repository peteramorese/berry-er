#pragma once

#include "HyperRectangle.h"

#include "berry/Options.h"

#include <Eigen/Core>

template <std::size_t DIM>
Eigen::MatrixXd BRY::HyperRectangle<DIM>::transformationMatrix(bry_deg_t m) const {
    BRY::bry_deg_t m_monoms = pow(m + 1, DIM);

    Eigen::MatrixXd T(m_monoms, m_monoms);
    T.setZero();

    Eigen::Vector<bry_float_t, DIM> scale = scaleFromUnit();
    Eigen::Vector<bry_float_t, DIM> translation = translationFromUnit();

    for (auto row_midx = mIdxW(DIM, m + 1); !row_midx.last(); ++row_midx) {
        std::vector<bry_idx_t> index_bounds(row_midx.size());
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
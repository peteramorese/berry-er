#pragma once

#include "Transform.h"

#include "berry/MultiIndex.h"
#include "berry/Operations.h"

template <std::size_t DIM>
Eigen::MatrixXd BRY::Transform::affine(const HyperRectangle<DIM>& set, bry_deg_t m) {
    BRY::bry_deg_t m_monoms = pow(m + 1, DIM);

    Eigen::MatrixXd T(m_monoms, m_monoms);
    T.setZero();

    Eigen::Vector<bry_float_t, DIM> scale = set.scaleFromUnit();
    Eigen::Vector<bry_float_t, DIM> translation = set.translationFromUnit();

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

            T(col_midx.inc().wrappedIdx(), row_midx.inc().wrappedIdx() + col_midx.inc().wrappedIdx()) = element;
        }
    }
    return T;
}
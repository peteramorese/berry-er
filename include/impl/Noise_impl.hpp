#pragma once

#include "Noise.h"

#include "berry/MultiIndex.h"
#include "berry/Operations.h"

template <std::size_t DIM>
BRY::Additive2ndMomentNoise<DIM>::Additive2ndMomentNoise(const Covariance<DIM>& covariance) 
    : m_cov(covariance)
{}

template <std::size_t DIM>
Eigen::MatrixXd BRY::Additive2ndMomentNoise<DIM>::additiveNoiseMatrix(bry_deg_t m) const {
    BRY::bry_deg_t m_monoms = pow(m + 1, DIM);

    Eigen::MatrixXd Gamma(m_monoms, m_monoms);
    Gamma.setZero();

    for (auto col_midx = mIdxW(DIM, m + 1); !col_midx.last(); ++col_midx) {

        std::vector<bry_idx_t> index_bounds(col_midx.size());
        for (std::size_t d = 0; d < DIM; ++d)
            index_bounds[d] = col_midx[d] + 1;

        for (auto row_midx = mIdxBEW(index_bounds, m + 1); !row_midx.last(); ++row_midx) {

            bry_float_t element = 1.0;

            bool zero = false;
            bry_deg_t moment = 0;
            int64_t first_cov_idx = -1;
            for (std::size_t j = 0; j < DIM; ++j) {
                bry_deg_t exp = col_midx[j] - row_midx[j];
                if (moment + exp > 2) {
                    zero = true;
                    break;
                } else if (exp == 2) {
                    element *= m_cov(j, j);
                    moment = 2;
                } else if (exp == 1) {
                    if (first_cov_idx == -1) {
                        first_cov_idx = j;
                        moment = 1;
                    } else {
                        element *= m_cov(first_cov_idx, j);
                        moment = 2;
                    }
                }
                element *= binom(col_midx[j], row_midx[j]);
            }

            if (!(zero || moment == 1)) {
                Gamma(row_midx.inc().wrappedIdx(), col_midx.inc().wrappedIdx()) = element;
            }
        }
    }
    
    return Gamma;
}

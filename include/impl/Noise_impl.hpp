#pragma once

#include "Noise.h"

#include "berry/MultiIndex.h"
#include "berry/Operations.h"

template <std::size_t DIM>
BRY::AdditiveGaussianNoise<DIM>::AdditiveGaussianNoise(const Covariance<DIM>& covariance) 
    : m_mean(Eigen::Vector<bry_float_t, DIM>::Zero())
    , m_cov(covariance)
{}

template <std::size_t DIM>
BRY::AdditiveGaussianNoise<DIM>::AdditiveGaussianNoise(const Eigen::Vector<bry_float_t, DIM>& mean, const Covariance<DIM>& covariance) 
    : m_mean(mean)
    , m_cov(covariance)
{}

template <std::size_t DIM>
BRY::Matrix BRY::AdditiveGaussianNoise<DIM>::additiveNoiseMatrix(bry_int_t m) const {
    BRY::bry_int_t m_monoms = pow(m + 1, DIM);

    Matrix Gamma(m_monoms, m_monoms);
    Gamma.setZero();

    for (auto col_midx = mIdxW(DIM, m + 1); !col_midx.last(); ++col_midx) {

        std::vector<bry_int_t> index_bounds(col_midx.size());
        for (std::size_t d = 0; d < DIM; ++d)
            index_bounds[d] = col_midx[d] + 1;

        for (auto row_midx = mIdxBEW(index_bounds, m + 1); !row_midx.last(); ++row_midx) {

            bry_float_t element = 1.0;

            bool zero = false;
            bry_int_t moment = 0;
            int64_t first_cov_idx = -1;
            for (std::size_t j = 0; j < DIM; ++j) {
                bry_int_t exp = col_midx[j] - row_midx[j];
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

template <std::size_t DIM>
BRY::Matrix BRY::AdditiveGaussianNoise<DIM>::momentMatrix(bry_int_t m) const {
    Polynomial<DIM> p(2);

    std::array<bry_int_t, DIM> mu_inds = makeUniformArray<bry_int_t, DIM>(0);
    for (std::size_t i = 0; i < DIM; ++i) {
        // Set the linear terms using the mean
        mu_inds[i] = 1;
        p.coeff(mu_inds) = m_mean[i];
        mu_inds[i] = 0;

        for (std::size_t j = 0; j < DIM; ++j) {
            std::array<bry_int_t, DIM> sigma_inds = makeUniformArray<bry_int_t, DIM>(0);
            sigma_inds[i] += 1;
            sigma_inds[j] += 1;
            p.coeff(sigma_inds) += 0.5 * m_cov(i, j);
        }
    }

    auto []<std::size_t SUB_DIM>(const Polynomial<DIM>& corner_polynomial, const std::array<bry_int_t, SUB_DIM>& span_dims) {
        std::vector<Polynomial<DIM>> cache_polynomials;
        cache_polynomials.reserve(SUB_DIM);
        for (std::size_t i = 0; i < SUB_DIM; i) {
            cache_polynomials.emplace_back(p);
        }
    }

    std::vector<Polynomial<DIM>> cache_polynomials;
    cache_polynomials.reserve(DIM);
    for (std::size_t i = 0; i < DIM; i) {
        cache_polynomials.emplace_back(p);
    }

    return Matrix(2,2);
}
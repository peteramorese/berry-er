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

    MomentGenerator gen(this, m);

    return Matrix(2,2);
}

template <std::size_t DIM>
BRY::AdditiveGaussianNoise<DIM>::MomentGenerator::MomentGenerator(const AdditiveGaussianNoise* enclosing, bry_int_t m) 
    : m_enclosing(enclosing)
    , m_m(m)
    , m_tensor(makeUniformArray<bry_int_t, DIM>(m + 1))
{
    Polynomial<DIM> p_exp(2);

    std::array<bry_int_t, DIM> mu_inds = makeUniformArray<bry_int_t, DIM>(0);
    for (std::size_t i = 0; i < DIM; ++i) {
        // Set the linear terms using the mean
        mu_inds[i] = 1;
        p_exp.coeff(mu_inds) = m_enclosing->m_mean[i];
        mu_inds[i] = 0;

        for (std::size_t j = 0; j < DIM; ++j) {
            std::array<bry_int_t, DIM> sigma_inds = makeUniformArray<bry_int_t, DIM>(0);
            sigma_inds[i] += 1;
            sigma_inds[j] += 1;
            p_exp.coeff(sigma_inds) += 0.5 * m_enclosing->m_cov(i, j);
        }
    }
    //Matrix(2,2)
    Eigen::Tensor<bry_float_t, DIM> moments(makeUniformArray<bry_int_t, DIM>(m + 1));

    std::array<bry_int_t, DIM> span_dims;
    for (std::size_t i = 0; i < DIM; ++i)
        span_dims[i] = i;

    DEBUG("ORIG Span dims ");
    for (auto e : span_dims)
        DEBUG(e);
    computeMoments<DIM>(p_exp, makeUniformArray<bry_int_t, DIM>(0), span_dims);
    
    bry_int_t tru_count = 0;
    for (auto midx = mIdx(DIM, m+1); !midx.last(); ++midx)
        ++tru_count;
    INFO("Count: " << m_counter << ", true count: " << tru_count);

}

template <std::size_t DIM>
template <std::size_t SUB_DIM>
void BRY::AdditiveGaussianNoise<DIM>::MomentGenerator::computeMoments(const Polynomial<DIM>& corner_coeff_polynomial, std::array<bry_int_t, DIM> corner_idx, const std::array<bry_int_t, SUB_DIM>& span_dims) { 
    std::string space(2*(DIM - SUB_DIM), ' ');
    //std::string space(1, ' ');
    std::cout << space << "SUBDIM: " << SUB_DIM << std::endl;
    std::cout << space << "span dims: [";
    for (std::size_t i = 0; i < SUB_DIM; ++i) {
        std::cout << span_dims[i] << ", ";
    }
    std::cout << "]" << std::endl;
    std::cout << space << "corner idx: ";
    for (std::size_t i = 0; i < DIM; ++i) {
        std::cout << corner_idx[i] << ", ";
    }
    std::cout << std::endl;

    if constexpr (SUB_DIM == 1) {
        ++corner_idx[span_dims[0]];
        for (; corner_idx[span_dims[0]] < m_m + 1; ++corner_idx[span_dims[0]]) {
            std::cout << "n ";
            for (std::size_t i = 0; i < DIM; ++i) {
                std::cout << corner_idx[i] << ", ";
            }
            std::cout << std::endl;
            ++m_counter;
        }
    } else {
            //DEBUG("NEW Span dims ");
            //for (auto e : new_span_dims)
            //    DEBUG(e);

        for (std::size_t corner_i = corner_idx[span_dims[0]]; corner_i < m_m + 1; ++corner_i) {
            
            // Copy the corner polynomial in to a cache polynomial for each subdim
            std::vector<Polynomial<DIM>> cache_polynomials;
            cache_polynomials.reserve(SUB_DIM);
            for (std::size_t i = 0; i < SUB_DIM; ++i) {
                cache_polynomials.emplace_back(corner_coeff_polynomial);
            }

            // Increment the corner indices along the span dimensions
            std::array<bry_int_t, DIM> new_corner_idx = corner_idx;
            for (bry_int_t span_idx : span_dims) {
                new_corner_idx[span_idx] = corner_i;
            }

            // If not the first corner after a split, then compute the coefficient
            if (corner_i > corner_idx[span_dims[0]]) {
                std::cout << "c ";
                for (std::size_t d = 0; d < DIM; ++d) {
                    std::cout << new_corner_idx[d] << ", ";
                }
                std::cout << std::endl;
                ++m_counter;
            }

            // When it makes it to the last index, don't bother making any recursive calls
            if (corner_i == m_m)
                continue;

            // For each remaining dimension, compute the moments for the perpendicular subspace
            for (std::size_t perp_dim_i = 0; perp_dim_i < SUB_DIM; ++perp_dim_i) {
                std::array<bry_int_t, SUB_DIM - 1> new_span_dims;
                auto it = new_span_dims.begin();
                for (std::size_t sub_dim_i = 0; sub_dim_i < SUB_DIM; ++sub_dim_i) {
                    if (sub_dim_i != perp_dim_i)
                        *(it++) = span_dims[sub_dim_i];
                }

                computeMoments<SUB_DIM - 1>(cache_polynomials[perp_dim_i], new_corner_idx, new_span_dims);

            }
        }
    }
};
#pragma once

#include "Noise.h"

#include "berry/MultiIndex.h"
#include "berry/Operations.h"

#include <algorithm>

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
    , m_duplicate(m_tensor.size(), false)
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

    struct ParDerivIncr {
        ParDerivIncr(const std::array<bry_int_t, DIM>& init_arr) : arr(init_arr) {}
        ParDerivIncr(const ParDerivIncr&) = default;
        void increment(bry_int_t incr_idx) {
            ++arr[incr_idx];
        }
        void compute(const std::array<bry_int_t, DIM>& idx) {
            std::cout << ">>";
            for (std::size_t i = 0; i < DIM; ++i) {
                std::cout << arr[i] << ", ";
            }
            std::cout << "\n";
        }
        std::array<bry_int_t, DIM> arr;
    };
    ParDerivIncr init_incr_obj(makeUniformArray<bry_int_t, DIM>(0));

    computeMoments<ParDerivIncr, DIM>(&init_incr_obj, makeUniformArray<bry_int_t, DIM>(0), span_dims);
    
    bry_int_t tru_count = 0;
    for (auto midx = mIdx(DIM, m+1); !midx.last(); ++midx)
        ++tru_count;
    INFO("Count: " << m_counter << ", true count: " << tru_count);

}

template <std::size_t DIM>
template <typename INCR_OBJECT, std::size_t SUB_DIM>
void BRY::AdditiveGaussianNoise<DIM>::MomentGenerator::computeMoments(INCR_OBJECT* corner_incr_obj, std::array<bry_int_t, DIM> corner_idx, const std::array<bry_int_t, SUB_DIM>& span_dims) { 
    std::string space(2*(DIM - SUB_DIM), ' ');
    //PAUSE;
    //std::cout << space << "SUBDIM: " << SUB_DIM << std::endl;
    //std::cout << space << "span dims: [";
    //for (std::size_t i = 0; i < SUB_DIM; ++i) {
    //    std::cout << span_dims[i] << ", ";
    //}
    //std::cout << "]" << std::endl;
    //std::cout << space << "corner idx: ";
    //for (std::size_t i = 0; i < DIM; ++i) {
    //    std::cout << corner_idx[i] << ", ";
    //}
    //std::cout << std::endl;

    if constexpr (SUB_DIM == 1) {
        ++corner_idx[span_dims[0]];
        for (; corner_idx[span_dims[0]] < m_m + 1; ++corner_idx[span_dims[0]]) {
            if (hasNotBeenSeen(corner_idx)) {
                std::cout << "n ";
                for (std::size_t i = 0; i < DIM; ++i) {
                    std::cout << corner_idx[i] << ", ";
                }
                std::cout << std::endl;
                corner_incr_obj->increment(span_dims[0]);
                corner_incr_obj->compute(corner_idx);
                ++m_counter;
            }
        }
    } else {
        bry_int_t max_starting_idx = *std::max_element(corner_idx.begin(), corner_idx.end());
        for (std::size_t corner_i = 0; corner_i < m_m + 1 - max_starting_idx; ++corner_i) {
            //std::cout << space << "corner " << corner_i << "/" << m_m - corner_idx[span_dims[0]] << std::endl;
            
            // Increment the corner indices along the span dimensions
            std::array<bry_int_t, DIM> new_corner_idx = corner_idx;
            for (bry_int_t span_idx : span_dims) {
                new_corner_idx[span_idx] += corner_i;
                if (corner_i)
                    corner_incr_obj->increment(span_idx);
            }

            // Copy the incr object into new caches
            std::vector<INCR_OBJECT> cache_incr_objs;
            cache_incr_objs.reserve(SUB_DIM);
            for (std::size_t i = 0; i < SUB_DIM; ++i) {
                cache_incr_objs.emplace_back(*corner_incr_obj);
            }

            //std::cout << "FIRST new corner idx 0: " << new_corner_idx[0] << ", 1: "<< new_corner_idx[1] << ", 2: "<< new_corner_idx[2] << std::endl;

            // If not the first corner after a split, then compute the coefficient
            //if (corner_i > 0) {
            if (hasNotBeenSeen(new_corner_idx)) {
                std::cout << "c ";
                for (std::size_t d = 0; d < DIM; ++d) {
                    std::cout << new_corner_idx[d] << ", ";
                }
                std::cout << std::endl;

                corner_incr_obj->compute(new_corner_idx);
                ++m_counter;
            }
            //}

            //// When it makes it to the last index, don't bother making any recursive calls
            if (corner_i == m_m - max_starting_idx)
                continue;

            // For each remaining dimension, compute the moments for the perpendicular subspace
            for (bry_int_t subspace_dim_i = 0; subspace_dim_i < SUB_DIM; ++subspace_dim_i) {

                // Offset the corner to prevent duplicates
                //bry_int_t inc_idx;
                //if constexpr (SUB_DIM > 2) {
                //    inc_idx = (subspace_dim_i) % SUB_DIM;
                //    //std::cout << space << "incrementing corner along idx " << inc_idx << std::endl;
                //    //PRINT_VEC3("new corner before inc ", new_corner_idx);
                //    ++new_corner_idx[inc_idx];
                //}

                //std::cout << "B4 SECOND new corner idx 0: " << new_corner_idx[0] << ", 1: "<< new_corner_idx[1] << ", 2: "<< new_corner_idx[2] << std::endl;
                std::array<bry_int_t, SUB_DIM - 1> new_span_dims;
                auto it = new_span_dims.begin();
                for (bry_int_t sub_dim_i = 0; sub_dim_i < SUB_DIM - 1; ++sub_dim_i) {
                    new_span_dims[sub_dim_i] = span_dims[(subspace_dim_i + sub_dim_i) % SUB_DIM];
                }

                //PRINT_VEC3("new corner before rec call ", new_corner_idx);
                computeMoments<INCR_OBJECT, SUB_DIM - 1>(&cache_incr_objs[subspace_dim_i], new_corner_idx, new_span_dims);

                // Remove the offset to make the index return to normal
                //if constexpr (SUB_DIM > 2) {
                //    --new_corner_idx[inc_idx];
                //}
            }
        }
    }
}

template <std::size_t DIM>
bool BRY::AdditiveGaussianNoise<DIM>::MomentGenerator::hasNotBeenSeen(const std::array<bry_int_t, DIM>& idx) {
    bry_int_t flat_idx = &m_tensor(idx) - m_tensor.data();
    if (m_duplicate[flat_idx]) {
        return false;
    } else {
        m_duplicate[flat_idx] = true;
        return true;
    }
}
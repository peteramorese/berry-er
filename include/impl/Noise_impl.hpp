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

    Eigen::Tensor<bry_float_t, DIM> moment_tensor = momentTensor(m);

    for (auto col_midx = mIdxW(DIM, m + 1); !col_midx.last(); ++col_midx) {

        std::vector<bry_int_t> index_bounds(col_midx.size());
        for (std::size_t d = 0; d < DIM; ++d)
            index_bounds[d] = col_midx[d] + 1;

        for (auto row_midx = mIdxBEW(index_bounds, m + 1); !row_midx.last(); ++row_midx) {

            bry_float_t element = 1.0;

            std::array<bry_int_t, DIM> moment_idx;
            bry_int_t moment_exp_sum = 0;
            for (bry_int_t j = 0; j < DIM; ++j) {
                moment_idx[j] = col_midx[j] - row_midx[j];
                moment_exp_sum += moment_idx[j];
                element *= binom(col_midx[j], row_midx[j]);
            }

            // Force remove odd-power moments for better numerical accuracy
            if (moment_exp_sum % 2 == 0) {
                element *= moment_tensor(moment_idx);
                Gamma(row_midx.inc().wrappedIdx(), col_midx.inc().wrappedIdx()) = element;
            }
        }
    }
    
    return Gamma.transpose();
}

template <std::size_t DIM>
Eigen::Tensor<BRY::bry_float_t, DIM> BRY::AdditiveGaussianNoise<DIM>::momentTensor(bry_int_t m) const {
    MomentGenerator gen(this, m);
    return gen.getTensor();
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
    //DEBUG("p_exp: " << p_exp);

    std::array<bry_int_t, DIM> span_dims;
    for (std::size_t i = 0; i < DIM; ++i)
        span_dims[i] = i;

    // Compute the first derivatives once and store them to save time
    std::vector<Polynomial<DIM>> exp_polynomial_first_derivatives;
    exp_polynomial_first_derivatives.reserve(DIM);
    for (bry_int_t d = 0; d < DIM; ++d) {
        exp_polynomial_first_derivatives.emplace_back(p_exp.derivative(d));
        //DEBUG("   p_exp deriv wrt x" << d << ": " << exp_polynomial_first_derivatives.back());
    }

    // Define the object to increment
    struct ParDerivIncr {
        ParDerivIncr(Eigen::Tensor<bry_float_t, DIM>* tensor_ptr_, 
                    const Polynomial<DIM>* exp_polynomial_, 
                    const std::vector<Polynomial<DIM>>* exp_polynomial_derivatives_, 
                    const Polynomial<DIM>& init_polynomial) 
            : tensor_ptr(tensor_ptr_) 
            , exp_polynomial(exp_polynomial_)
            , exp_polynomial_derivatives(exp_polynomial_derivatives_)
            , polynomial(init_polynomial)
            {}
        ParDerivIncr(const ParDerivIncr&) = default;

        void increment(bry_int_t incr_idx) {
            // Take the partial derivative of the MGF wrt the incr idx

            //DEBUG("curr poly: " << polynomial << " deg: " << polynomial.degree());
            //DEBUG("curr poly derivative: " << polynomial.derivative(incr_idx));
            //DEBUG("exp poly derivative:  " << (*exp_polynomial_derivatives)[incr_idx] << " deg: " << (*exp_polynomial_derivatives)[incr_idx].degree());
            //DEBUG("prod:  " << polynomial * (*exp_polynomial_derivatives)[incr_idx]);
            polynomial = polynomial.derivative(incr_idx) + polynomial * (*exp_polynomial_derivatives)[incr_idx];
        }

        void compute(const std::array<bry_int_t, DIM>& idx) {
            (*tensor_ptr)(idx) = (*polynomial.tensor().data());
        }

        Eigen::Tensor<bry_float_t, DIM>* tensor_ptr;
        const Polynomial<DIM>* exp_polynomial;
        const std::vector<Polynomial<DIM>>* exp_polynomial_derivatives;
        Polynomial<DIM> polynomial;
    };

    // Make the initial object
    Polynomial<DIM> init_polynomial(0);
    init_polynomial.coeff(makeUniformArray<bry_int_t, DIM>(0)) = 1.0;
    //DEBUG("Init polynomial: " << init_polynomial);
    ParDerivIncr init_incr_obj(&m_tensor, &p_exp, &exp_polynomial_first_derivatives, init_polynomial);

    INFO("Calculating moment matrix...");
    minPathIncr<ParDerivIncr, DIM>(&init_incr_obj, makeUniformArray<bry_int_t, DIM>(0), span_dims);
    NEW_LINE;
    INFO("Done!");
}

template <std::size_t DIM>
const Eigen::Tensor<BRY::bry_float_t, DIM>& BRY::AdditiveGaussianNoise<DIM>::MomentGenerator::getTensor() const {
    return m_tensor;
}

template <std::size_t DIM>
template <typename INCR_OBJECT, std::size_t SUB_DIM>
void BRY::AdditiveGaussianNoise<DIM>::MomentGenerator::minPathIncr(INCR_OBJECT* corner_incr_obj, std::array<bry_int_t, DIM> corner_idx, const std::array<bry_int_t, SUB_DIM>& span_dims) { 
    //std::string space(2*(DIM - SUB_DIM), ' ');
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
                //std::cout << "n ";
                //for (std::size_t i = 0; i < DIM; ++i) {
                //    std::cout << corner_idx[i] << ", ";
                //}
                //std::cout << std::endl;
                corner_incr_obj->increment(span_dims[0]);
                corner_incr_obj->compute(corner_idx);
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
                //std::cout << "c ";
                //for (std::size_t d = 0; d < DIM; ++d) {
                //    std::cout << new_corner_idx[d] << ", ";
                //}
                //std::cout << std::endl;

                corner_incr_obj->compute(new_corner_idx);
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
                minPathIncr<INCR_OBJECT, SUB_DIM - 1>(&cache_incr_objs[subspace_dim_i], new_corner_idx, new_span_dims);

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
        ++m_counter;
        INFO_SMLN(" " << m_counter << " / " << m_tensor.size());
        m_duplicate[flat_idx] = true;
        return true;
    }
}
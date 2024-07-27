#pragma once

#include "Dynamics.h"

#include "berry/MultiIndex.h"

#include "lemon/Random.h"

#include <iomanip>

template <std::size_t DIM>
template <typename ... DEGS>
BRY::PolynomialDynamics<DIM>::PolynomialDynamics(DEGS ... degrees) 
{
    static_assert(is_uniform_convertible_type<bry_int_t, DEGS...>(), "All parameters passed to ctor must be convertible to bry_int_t");
    ExponentVec<DIM> degree_vec = makeExponentVec(degrees...);
    m_f.reserve(DIM);
    for (bry_int_t deg : degree_vec)
        m_f.emplace_back(deg);
}

template <std::size_t DIM>
BRY::PolynomialDynamics<DIM>::PolynomialDynamics(std::array<bry_int_t, DIM> degrees)
{
    m_f.reserve(DIM);
    for (bry_int_t deg : degrees)
        m_f.emplace_back(deg);
}

template <std::size_t DIM>
BRY::Polynomial<DIM, BRY::Basis::Power>& BRY::PolynomialDynamics<DIM>::operator[](std::size_t j) {
    return m_f[j];
}

template <std::size_t DIM>
const BRY::Polynomial<DIM, BRY::Basis::Power>& BRY::PolynomialDynamics<DIM>::operator[](std::size_t j) const {
    return m_f[j];
}

template <std::size_t DIM>
std::array<BRY::bry_int_t, DIM> BRY::PolynomialDynamics<DIM>::degrees() const {
    std::array<bry_int_t, DIM> deg_array;
    for (std::size_t i = 0; i < DIM; ++i) {
        deg_array[i] = m_f[i].degree();
    }
    return deg_array;
}

template <std::size_t DIM>
BRY::bry_int_t BRY::PolynomialDynamics<DIM>::summedDegree() const {
    bry_int_t n_sum = 0;
    for (const auto& polynomial : m_f)
        n_sum += polynomial.degree();
    return n_sum;
}

template <std::size_t DIM>
BRY::bry_int_t BRY::PolynomialDynamics<DIM>::composedDegree(bry_int_t m) const {
    return m * summedDegree();
}

template <std::size_t DIM>
void BRY::PolynomialDynamics<DIM>::setRandom() {
    for (bry_int_t i = 0; i < DIM; ++i) {
        std::array<bry_int_t, DIM> idx;
        MultiIndex<ExhaustiveIncrementer> midx(idx.data(), DIM, true, m_f[i].degree() + 1);
        for (; !midx.last(); ++midx) {
            bool diag_lin = true;
            for (bry_int_t j = 0; j < DIM; ++j) {
                if ((j == i && midx[j] != 1) || (j != i && midx[j] != 0)) {
                    diag_lin = false;
                    break;
                } 
            }
            if (diag_lin) {
                m_f[i].coeff(idx) = lemon::RNG::srandd(0.9, 1.1);
            } else {
                m_f[i].coeff(idx) = lemon::RNG::srandd(-0.1, 0.1);
            }
        }
    }
}

template <std::size_t DIM>
BRY::Matrix BRY::PolynomialDynamics<DIM>::dynamicsPowerMatrix(bry_int_t m) const {
    BRY::bry_int_t p = composedDegree(m);
    BRY::bry_int_t m_monoms = BRY::pow(m + 1, DIM);
    BRY::bry_int_t p_monoms = BRY::pow(p + 1, DIM);

    Matrix F(p_monoms, m_monoms);
    F.col(0) = Vector::Zero(p_monoms);
    F(0, 0) = 1.0;

    std::array<bry_int_t, DIM> dimensions;
    for (std::size_t d = 0; d < DIM; ++d)
        dimensions[d] = d;

    std::array<Eigen::Tensor<bry_complex_t, DIM>, DIM> complex_tensors;
    for (std::size_t d = 0; d < DIM; ++d) {
        Eigen::Tensor<bry_float_t, DIM> padded_tensor = _BRY::expandToMatchSize<DIM>(m_f[d].tensor(), p + 1);;
        complex_tensors[d] = padded_tensor.template fft<Eigen::BothParts, Eigen::FFT_FORWARD>(dimensions);
    }

    // Temporary variable to hold the product of each dynamics polynomial
    Eigen::Tensor<bry_complex_t, DIM> product(makeUniformArray<bry_int_t, DIM>(p + 1));

    auto col_midx = mIdxW(DIM, m + 1);
    ++col_midx;
    for (; !col_midx.last(); ++col_midx) {
        product.setConstant(1.0);
        for (std::size_t d = 0; d < DIM; ++d) {
            bry_int_t exponent = col_midx[d];
            if (exponent > 0)
                product *= complex_tensors[d].pow(exponent);
        }

        Eigen::TensorMap<Eigen::Tensor<double, 1>> column(F.col(col_midx.inc().wrappedIdx()).data(), p_monoms);

        Eigen::Tensor<BRY::bry_float_t, DIM> product_coefficients = product.template fft<Eigen::RealPart, Eigen::FFT_REVERSE>(dimensions);

        std::array<bry_int_t, 1> one_dim{{p_monoms}};
        column = product_coefficients.reshape(one_dim);
    }

    return F;
}

template <std::size_t DIM>
std::ostream& operator<<(std::ostream& os, const BRY::PolynomialDynamics<DIM>& dynamics) {
    BRY::bry_int_t max_header_sz = 0;
    std::array<std::string, DIM> headers;
    for (BRY::bry_int_t i = 0; i < DIM; ++i) {
        headers[i] = "x" + std::to_string(i) + "'";
        if (headers[i].size() > max_header_sz) {
            max_header_sz = headers[i].size();
        }
    }
    for (BRY::bry_int_t i = 0; i < DIM; ++i) {
        os << std::setw(max_header_sz) << LMN_LOG_BWHITE(headers[i]) << " = " << dynamics[i] << "\n";
    }
    return os;
}
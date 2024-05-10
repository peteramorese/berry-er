#pragma once

#include "Dynamics.h"

#include "berry/MultiIndex.h"

template <std::size_t DIM>
template <typename ... DEGS>
BRY::PolynomialDynamics<DIM>::PolynomialDynamics(DEGS ... degrees) 
{
    static_assert(is_uniform_convertible_type<bry_deg_t, DEGS...>(), "All parameters passed to ctor must be convertible to bry_deg_t");
    ExponentVec<DIM> degree_vec = makeExponentVec(degrees...);
    m_f.reserve(DIM);
    for (bry_deg_t deg : degree_vec)
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
std::array<BRY::bry_deg_t, DIM> BRY::PolynomialDynamics<DIM>::degrees() const {
    std::array<bry_deg_t, DIM> deg_array;
    for (std::size_t i = 0; i < DIM; ++i) {
        deg_array[i] = m_f[i].degree();
    }
    return deg_array;
}

template <std::size_t DIM>
BRY::bry_deg_t BRY::PolynomialDynamics<DIM>::summedDegree() const {
    bry_deg_t n_sum = 0;
    for (const auto& polynomial : m_f)
        n_sum += polynomial.degree();
    return n_sum;
}

template <std::size_t DIM>
BRY::bry_deg_t BRY::PolynomialDynamics<DIM>::composedDegree(bry_deg_t m) const {
    return m * summedDegree();
}

template <std::size_t DIM>
Eigen::MatrixXd BRY::PolynomialDynamics<DIM>::dynamicsPowerMatrix(bry_deg_t m) const {
    BRY::bry_deg_t p = composedDegree(m);
    BRY::bry_deg_t m_monoms = pow(m + 1, DIM);
    BRY::bry_deg_t p_monoms = pow(p + 1, DIM);

    Eigen::MatrixXd F(p_monoms, m_monoms);
    F.col(0) = Eigen::VectorXd::Zero(p_monoms);
    F(0, 0) = 1.0;

    std::array<bry_idx_t, DIM> dimensions;
    for (std::size_t d = 0; d < DIM; ++d)
        dimensions[d] = d;

    std::array<Eigen::Tensor<bry_complex_t, DIM>, DIM> complex_tensors;
    for (std::size_t d = 0; d < DIM; ++d) {
        Eigen::Tensor<bry_float_t, DIM> padded_tensor = _BRY::expandToMatchSize<DIM>(m_f[d].tensor(), p + 1);;
        complex_tensors[d] = padded_tensor.template fft<Eigen::BothParts, Eigen::FFT_FORWARD>(dimensions);
    }

    // Temporary variable to hold the product of each dynamics polynomial
    Eigen::Tensor<bry_complex_t, DIM> product(makeUniformArray<bry_deg_t, DIM>(p + 1));

    auto col_midx = mIdxW(DIM, m + 1);
    ++col_midx;
    for (; !col_midx.last(); ++col_midx) {
        product.setConstant(1.0);
        for (std::size_t d = 0; d < DIM; ++d) {
            bry_deg_t exponent = col_midx[d];
            if (exponent > 0)
                product *= complex_tensors[d].pow(exponent);
        }

        Eigen::TensorMap<Eigen::Tensor<double, 1>> column(F.col(col_midx.inc().wrappedIdx()).data(), p_monoms);

        Eigen::Tensor<BRY::bry_float_t, DIM> product_coefficients = product.template fft<Eigen::RealPart, Eigen::FFT_REVERSE>(dimensions);

        std::array<bry_deg_t, 1> one_dim{{p_monoms}};
        column = product_coefficients.reshape(one_dim);
    }

    return F;
}
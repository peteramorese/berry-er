#pragma once

#include "Dynamics.h"

template <std::size_t DIM>
BRY::Polynomial<DIM, Basis::Power>& BRY::PolynomialDynamics::operator[](std::size_t j) {
    return m_f[j];
}

template <std::size_t DIM>
const BRY::Polynomial<DIM, Basis::Power>& BRY::PolynomialDynamics::operator[](std::size_t j) const {
    return m_f[j];
}

template <std::size_t DIM>
BRY_INL std::array<bry_deg_t, DIM> BRY::PolynomialDynamics::degrees() const {
    std::array<bry_deg_t, DIM> deg_array;
    for (std::size_t i = 0; i < DIM; ++i) {
        deg_array[i] = m_f[i].degree();
    }
    return deg_array;
}

template <std::size_t DIM>
Eigen::MatrixXd BRY::PolynomialDynamics::dynamicsPowerMatrix(bry_deg_t m) const {
    auto deg_arr = degrees();
    BRY::bry_deg_t n_sum = std::accumulate(deg_arr.begin(), deg_arr.end(), 0);
    BRY::bry_deg_t p = m * n_sum;
    BRY::bry_deg_t m_monoms = pow(m, DIM) + 1;
    BRY::bry_deg_t p_monoms = pow(p, DIM) + 1;

    Eigen::MatrixXd F(p_monoms, m_monoms);
    F.col(0) = Eigen::VectorXd::Ones(p_monoms);

    std::array<bry_idx_t, DIM> dimensions;
    for (std::size_t d = 0; d < DIM; ++d)
        dimensions[d] = d;

    std::array<Eigen::Tensor<bry_complex_t, DIM>, DIM> complex_tensors;
    for (std::size_t d = 0; d < DIM; ++d) {
        Eigen::Tensor<bry_complex_t, DIM> padded_tensor = _BRY::expandToMatchSize<DIM>(m_f[d], p + 1);;
        complex_tensors[d] = padded_tensor.template fft<Eigen::BothParts, Eigen::FFT_FORWARD>(dimensions);
    }

    // Temporary variable to hold the product of each dynamics polynomial
    Eigen::Tensor<bry_complex_t, DIM> t_product(makeUniformArray<bry_deg_t, DIM>(p + 1));

    for (auto col_midx = mIdxW(DIM, m + 1); !col_midx.last(); ++col_midx) {
        product.setOnes();
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
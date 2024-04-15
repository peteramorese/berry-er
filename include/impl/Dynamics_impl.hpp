#pragma once

#include "Dynamics.h"

template <std::size_t DIM>

BRY::Polynomial<DIM, Basis::Power>& BRY::PolynomialDynamics::operator[](std::size_t j) {
    return m_f[j];
}

const BRY::Polynomial<DIM, Basis::Power>& BRY::PolynomialDynamics::operator[](std::size_t j) const {
    return m_f[j];
}
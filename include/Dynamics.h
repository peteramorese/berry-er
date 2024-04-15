#pragma once

#include "berry/Polynomial.h"

#include <array>

namespace BRY {

//template <std::size_t DIM, Basis BASIS = Basis::Power>
//using PolyPtr = std::shared_ptr<Polynomial<DIM, BASIS>>;

template <std::size_t DIM>
class PolynomialDynamics {
    public:
        PolynomialDynamics() {}
    
        BRY_INL Polynomial<DIM, Basis::Power>& operator[](std::size_t j);
        BRY_INL const Polynomial<DIM, Basis::Power>& operator[](std::size_t j) const;

    private:
        std::array<Polynomial<DIM, Basis::Power>, DIM> m_f;

};

/* TODO */
class BoundedPolynomialDynamics {};

}

#include "impl/Dynamics_impl.hpp"
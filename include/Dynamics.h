#pragma once

#include "berry/Polynomial.h"
#include "berry/Types.h"

#include <array>

#include <Eigen/Dense>

namespace BRY {

template <std::size_t DIM>
class PolynomialDynamics {
    public:
        PolynomialDynamics() = delete;
        
        template <typename ... DEGS>
        PolynomialDynamics(DEGS ... degrees);
    
        BRY_INL Polynomial<DIM, Basis::Power>& operator[](std::size_t j);
        BRY_INL const Polynomial<DIM, Basis::Power>& operator[](std::size_t j) const;

        BRY_INL std::array<bry_int_t, DIM> degrees() const;
        BRY_INL bry_int_t summedDegree() const; 
        BRY_INL bry_int_t composedDegree(bry_int_t m) const; 

        /// @brief Compute the dynamics power matrix 'F'
        /// @param m Max exponent of the dynamics (degree of barrier)
        /// @return Dynamics composition matrix of size `p^d` by `m^d` 
        Matrix dynamicsPowerMatrix(bry_int_t m) const;

    private:
        std::vector<Polynomial<DIM, Basis::Power>> m_f;

};

/* TODO */
class BoundedPolynomialDynamics {};

}

#include "impl/Dynamics_impl.hpp"
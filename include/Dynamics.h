#pragma once

#include "berry/Polynomial.h"

#include <array>

namespace BRY {

//template <std::size_t DIM, Basis BASIS = Basis::Power>
//using PolyPtr = std::shared_ptr<Polynomial<DIM, BASIS>>;

template <std::size_t DIM>
class PolynomialDynamics {
    public:
        PolynomialDynamics() = delete;
        
        template <typename ... DEGS>
        PolynomialDynamics(DEGS ... degrees);
    
        BRY_INL Polynomial<DIM, Basis::Power>& operator[](std::size_t j);
        BRY_INL const Polynomial<DIM, Basis::Power>& operator[](std::size_t j) const;

        BRY_INL std::array<bry_deg_t, DIM> degrees() const;

        /// @brief Compute the dynamics power matrix 'F'
        /// @param m Max exponent of the dynamics (degree of barrier)
        /// @return Dynamics composition matrix of size `p^d` by `m^d` 
        Eigen::MatrixXd dynamicsPowerMatrix(bry_deg_t m) const;

    private:
        std::vector<Polynomial<DIM, Basis::Power>> m_f;

};

/* TODO */
class BoundedPolynomialDynamics {};

}

#include "impl/Dynamics_impl.hpp"
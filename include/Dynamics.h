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
        
        PolynomialDynamics(std::array<bry_int_t, DIM> degrees);
    
        BRY_INL Polynomial<DIM, Basis::Power>& operator[](std::size_t j);
        BRY_INL const Polynomial<DIM, Basis::Power>& operator[](std::size_t j) const;

        /// @brief Degrees of each polynomial component
        /// @return Array of degrees corresponding to each dimension
        BRY_INL std::array<bry_int_t, DIM> degrees() const;

        /// @brief Sum of all polynomial component degrees
        /// @return Summed degree
        BRY_INL bry_int_t summedDegree() const; 

        /// @brief Degree of the polynomial after being composed with a barrier of degree `m`
        /// @param m Barrier degree
        /// @return `m * summedDegree()`
        BRY_INL bry_int_t composedDegree(bry_int_t m) const; 

        /// @brief Make the dynamics random. Diagonal linear terms have coeffs btw [.9, 1.1], all others have coeffs btw [-.1, .1]
        void setRandom();

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

template <std::size_t DIM>
std::ostream& operator<<(std::ostream& os, const BRY::PolynomialDynamics<DIM>& dynamics);


#include "impl/Dynamics_impl.hpp"
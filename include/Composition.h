#pragma once

#include "berry/Polynomial.h"

#include "Dynamics.h"

#include <Eigen/Core>

namespace BRY {

/// @brief Composer object for computing the full composed matrix 'M' for polynomial dynamics
/// @tparam DIM 
template <std::size_t DIM>
class Composer {
    public:
        Composer(const std::shared_ptr<PolynomialDynamics<DIM>>& dynamics, bry_deg_t barrier_degree);

        /// @brief Compute the composition matrix 'M'
        /// @param covariance Covariance of the noise term 'v'
        /// @return Barrier composition matrix 'M'
        Eigen::MatrixXd compositionMatrix(const Eigen::Matrix<bry_float_t, DIM, DIM>& covariance);

    private:
        const std::shared_ptr<const PolynomialDynamics<DIM>> m_dynamics;
        bry_deg_t m_barrier_degree;
};

}

#include "impl/Composition_impl.h"
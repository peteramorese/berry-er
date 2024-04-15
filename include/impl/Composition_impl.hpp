#pragma once

#include "Composition.h"

template <std::size_t DIM>
BRY::Composer::Composer(const PolyDynamicsPtr<DIM>& dynamics, bry_deg_t barrier_degree) 
    : m_dynamics(dynamics)
    , m_barrier_degree(barrier_degree)
{}

Eigen::MatrixXd compositionMatrix(const Eigen::Matrix<bry_float_t, DIM, DIM>& covariance) {


}

        
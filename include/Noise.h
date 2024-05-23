#pragma once

#include <array>

#include <Eigen/Dense>

namespace BRY {

template <std::size_t DIM>
using Covariance = Eigen::Matrix<bry_float_t, DIM, DIM>;

template <std::size_t DIM>
class Additive2ndMomentNoise {
    public:
        Additive2ndMomentNoise() = delete;
        Additive2ndMomentNoise(const Covariance<DIM>& covariance);

        Eigen::MatrixXd additiveNoiseMatrix(bry_int_t m) const;
        
    private:
        Covariance<DIM> m_cov;

};

}

#include "impl/Noise_impl.hpp"
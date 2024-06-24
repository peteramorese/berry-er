#pragma once

#include <array>

#include <Eigen/Dense>

namespace BRY {

template <std::size_t DIM>
using Covariance = Eigen::Matrix<bry_float_t, DIM, DIM>;

template <std::size_t DIM>
class AdditiveGaussianNoise {
    public:
        AdditiveGaussianNoise() = delete;
        AdditiveGaussianNoise(const Covariance<DIM>& covariance);
        AdditiveGaussianNoise(const Eigen::Vector<bry_float_t, DIM>& mean, const Covariance<DIM>& covariance);

        Matrix additiveNoiseMatrix(bry_int_t m) const;

        Matrix momentMatrix(bry_int_t m) const;
        
    private:
        Eigen::Vector<bry_float_t, DIM> m_mean;
        Covariance<DIM> m_cov;

};

}

#include "impl/Noise_impl.hpp"
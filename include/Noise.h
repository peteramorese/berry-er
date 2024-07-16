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

        Eigen::Tensor<bry_float_t, DIM> momentTensor(bry_int_t m) const;

    private:
        class MomentGenerator {
            public:
                MomentGenerator(const AdditiveGaussianNoise* enclosing, bry_int_t m);

                const Eigen::Tensor<bry_float_t, DIM>& getTensor() const;

            private:
                template <typename INCR_OBJECT, std::size_t SUB_DIM>
                void minPathIncr(INCR_OBJECT* corner_incr_obj, std::array<bry_int_t, DIM> curr_idx, const std::array<bry_int_t, SUB_DIM>& span_dims);

                bool hasNotBeenSeen(const std::array<bry_int_t, DIM>& idx);

            private:
                const AdditiveGaussianNoise* m_enclosing;
                bry_int_t m_m;
                Eigen::Tensor<bry_float_t, DIM> m_tensor;
                bry_int_t m_counter = 0;
                std::vector<bool> m_duplicate;
        };
        
    private:
        Eigen::Vector<bry_float_t, DIM> m_mean;
        Covariance<DIM> m_cov;

};

}

#include "impl/Noise_impl.hpp"
#pragma once

#include "berry/Types.h"

#include <memory>

namespace BRY {


template <std::size_t DIM>
class MonomialFilter {
    public:
        /// @brief Determine whether or not to remove a monomial based on the exponent vector
        /// @param exponent_vec Exponent vector (contiguous array) of monomial in question. Size must be DIM
        /// @return `true` if filter should remove and `false` otherwise
        virtual bool remove(const bry_int_t* exponent_vec) = 0;
};

/// @brief Remove all elements with summed degree greater than the degree of the polynomial
template <std::size_t DIM>
class DiagDegFilter : public MonomialFilter<DIM> {
    public:
        DiagDegFilter(bry_int_t barrier_deg);
        virtual bool remove(const bry_int_t* exponent_vec) override;
    private:    
        bry_int_t m_barrier_deg;
};

/// @brief Remove all elements with odd summed degree
template <std::size_t DIM>
class OddSumFilter : public MonomialFilter<DIM> {
    public:
        virtual bool remove(const bry_int_t* exponent_vec) override;
};

/// @brief Combine filters such that the element is removed only if all filters remove the element
template <std::size_t DIM>
class ConjMultiFilter : public MonomialFilter<DIM> {
    public:
        ConjMultiFilter(const std::vector<std::shared_ptr<MonomialFilter<DIM>>>& filters);

        virtual bool remove(const bry_int_t* exponent_vec) override;
    private:
        std::vector<std::shared_ptr<MonomialFilter<DIM>>> m_filters;
};


}

#include "impl/MonomialFilter_impl.hpp"
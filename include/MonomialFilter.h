#pragma once

#include "berry/Types.h"

#include <memory>

namespace BRY {


template <std::size_t DIM>
class MonomialFilter {
    public:
        MonomialFilter(bry_int_t barrier_deg);

        /// @brief Determine whether or not to remove a monomial based on the exponent vector
        /// @param exponent_vec Exponent vector (contiguous array) of monomial in question. Size must be DIM
        /// @return `true` if filter should remove and `false` otherwise
        virtual bool remove(const bry_int_t* exponent_vec) const = 0;

        BRY_INL bry_int_t barrierDeg() const;

        /// @brief Get the number remaining monomials after the filter is applied
        /// @return Number of remaining monoms
        BRY_INL bry_int_t nRemainingMonoms() const; 

        /// @brief Flags indicating if a monomial should be removed or not
        /// @return Get a flattened array of flags where `true` corresponds to `remove(wrapped_idx) = true`
        BRY_INL const Eigen::Tensor<bool, DIM>& flags() const; 

        /// @brief Get the new wrapped index of the filtered polynomial given the old wrapped index of the unfiltered polynomial
        /// @param old_wrapped_idx Idx of unfiltered polynomial
        /// @return Idx of filtered polynomial
        BRY_INL const bry_int_t newWrappedIdx(bry_int_t old_wrapped_idx) const; 

        /// @brief Apply the filter and remove rows corresponding to filtered monomials
        /// @param matrix Matrix to be edited in-place
        Matrix applyToCoeffMatrixRows(const Matrix& matrix) const;

        /// @brief Apply the filter and remove cols corresponding to filtered monomials
        /// @param matrix Matrix to be edited in-place
        Matrix applyToCoeffMatrixCols(const Matrix& matrix) const;

    protected:
        /// @brief Must call this function in any derived class constructor
        void init();

    protected:
        Eigen::Tensor<bool, DIM> m_flags;
        std::vector<bry_int_t> m_idx_map;
        const bry_int_t m_barrier_deg; 
        bry_int_t m_remaining_monoms;
};

/// @brief Remove all elements with summed degree greater than the degree of the polynomial
template <std::size_t DIM>
class DiagDegFilter : public MonomialFilter<DIM> {
    public:
        DiagDegFilter(bry_int_t barrier_deg);
        virtual bool remove(const bry_int_t* exponent_vec) const override;
};

/// @brief Remove all elements with odd summed degree
template <std::size_t DIM>
class OddSumFilter : public MonomialFilter<DIM> {
    public:
        OddSumFilter(bry_int_t barrier_deg);
        OddSumFilter(bry_int_t barrier_deg, bry_int_t min_sum_exp_to_keep);
        virtual bool remove(const bry_int_t* exponent_vec) const override;

    private:
        bry_int_t m_min_sum_exp_to_keep = 3;
};

template <std::size_t DIM>
class UniformRandomFilter : public MonomialFilter<DIM> {
    public:
        UniformRandomFilter(bry_int_t barrier_deg, bry_int_t monoms_to_remove);
        virtual bool remove(const bry_int_t* exponent_vec) const override;

    private:
        bry_int_t m_monoms_to_remove;
};

///// @brief Combine filters such that the element is removed only if all filters remove the element
//template <std::size_t DIM>
//class ConjMultiFilter : public MonomialFilter<DIM> {
//    public:
//        ConjMultiFilter(const std::vector<std::shared_ptr<MonomialFilter<DIM>>>& filters);
//
//        virtual bool remove(const bry_int_t* exponent_vec) override;
//    private:
//        std::vector<std::shared_ptr<MonomialFilter<DIM>>> m_filters;
//};


}

#include "impl/MonomialFilter_impl.hpp"
#pragma once

#include "MonomialFilter.h"

#include "berry/MultiIndex.h"
#include "berry/Operations.h"
#include "lemon/Random.h"

#include <algorithm>

template <std::size_t DIM>
BRY::MonomialFilter<DIM>::MonomialFilter(bry_int_t barrier_deg) 
    : m_barrier_deg(barrier_deg)
    , m_flags(makeUniformArray<bry_int_t, DIM>(barrier_deg + 1))
    , m_idx_map(pow(barrier_deg + 1, DIM))
{}

template <std::size_t DIM>
void BRY::MonomialFilter<DIM>::init() {
    m_remaining_monoms = 0;
    bool* flags_arr = m_flags.data();
    DEBUG("flags size: " << m_flags.size());
    for (auto midx = mIdxW(DIM, m_barrier_deg + 1); !midx.last(); ++midx) {
        bool remove_monom = remove(midx.begin());
        flags_arr[midx.inc().wrappedIdx()] = remove_monom;
        m_idx_map[midx.inc().wrappedIdx()] = m_remaining_monoms;
        m_remaining_monoms += !remove_monom;
        DEBUG("i: " << midx.inc().wrappedIdx() << " remove monom: " << remove_monom);
    }
}

template <std::size_t DIM>
BRY::bry_int_t BRY::MonomialFilter<DIM>::barrierDeg() const {
    return m_barrier_deg;
}

template <std::size_t DIM>
BRY::bry_int_t BRY::MonomialFilter<DIM>::nRemainingMonoms() const {
    return m_remaining_monoms;
}

template <std::size_t DIM>
const Eigen::Tensor<bool, DIM>& BRY::MonomialFilter<DIM>::flags() const {
    return m_flags;
}

template <std::size_t DIM>
const BRY::bry_int_t BRY::MonomialFilter<DIM>::newWrappedIdx(bry_int_t old_wrapped_idx) const {
    return m_idx_map[old_wrapped_idx];
}

template <std::size_t DIM>
BRY::Matrix BRY::MonomialFilter<DIM>::applyToCoeffMatrixRows(const Matrix& matrix) const {
    ASSERT(matrix.rows() == m_flags.size(), "Number of rows in matrix does not match number monomials filter should consider");
    Matrix filtered_matrix(m_remaining_monoms, matrix.cols());
    bry_int_t new_idx = 0;
    for (bry_int_t old_idx = 0; old_idx < matrix.rows(); ++old_idx) {
        if (!m_flags.data()[old_idx]) {
            filtered_matrix.row(new_idx++) = matrix.row(old_idx);
        }
    }
    return filtered_matrix;
}

template <std::size_t DIM>
BRY::Matrix BRY::MonomialFilter<DIM>::applyToCoeffMatrixCols(const Matrix& matrix) const {
    ASSERT(matrix.cols() == m_flags.size(), "Number of cols in matrix does not match number monomials filter should consider");
    Matrix filtered_matrix(matrix.rows(), m_remaining_monoms);
    bry_int_t new_idx = 0;
    for (bry_int_t old_idx = 0; old_idx < matrix.cols(); ++old_idx) {
        if (!m_flags.data()[old_idx]) {
            filtered_matrix.col(new_idx++) = matrix.col(old_idx);
        }
    }
    return filtered_matrix;
}

template <std::size_t DIM>
BRY::DiagDegFilter<DIM>::DiagDegFilter(bry_int_t barrier_deg)
    : MonomialFilter<DIM>(barrier_deg)
{
    this->init();
}

template <std::size_t DIM>
bool BRY::DiagDegFilter<DIM>::remove(const bry_int_t* exponent_vec) const {
    return std::accumulate(exponent_vec, exponent_vec + DIM, 0) > this->m_barrier_deg;
}

template <std::size_t DIM>
BRY::OddSumFilter<DIM>::OddSumFilter(bry_int_t barrier_deg)
    : MonomialFilter<DIM>(barrier_deg)
{
    this->init();
}

template <std::size_t DIM>
BRY::OddSumFilter<DIM>::OddSumFilter(bry_int_t barrier_deg, bry_int_t min_sum_exp_to_keep)
    : MonomialFilter<DIM>(barrier_deg)
    , m_min_sum_exp_to_keep(min_sum_exp_to_keep)
{
    this->init();
}

template <std::size_t DIM>
bool BRY::OddSumFilter<DIM>::remove(const bry_int_t* exponent_vec) const {
    bry_int_t sum_exp = std::accumulate(exponent_vec, exponent_vec + DIM, 0);
    if (sum_exp < m_min_sum_exp_to_keep) {
        return false;
    }
    return (std::accumulate(exponent_vec, exponent_vec + DIM, 0) + 1) % 2 == 0;
}

template <std::size_t DIM>
BRY::UniformRandomFilter<DIM>::UniformRandomFilter(bry_int_t barrier_deg, bry_int_t monoms_to_remove)
    : MonomialFilter<DIM>(barrier_deg)
    , m_monoms_to_remove(monoms_to_remove)
{
    //this->init();
    std::vector<bool> seen(this->m_flags.size(), false);
    bry_int_t monoms_removed = 0;
    bool* flags_arr = this->m_flags.data();
    while (monoms_removed < m_monoms_to_remove) {
        bry_int_t idx_to_remove = lemon::RNG::randi((bry_int_t)1, (bry_int_t)this->m_flags.size());
        if (!seen[idx_to_remove]) {
            flags_arr[idx_to_remove] = true;
            seen[idx_to_remove] = true;
            ++monoms_removed;
        }
        DEBUG("hello");
    }
    INFO("Done creating random filter");

    this->init();
    DEBUG("remaining monoms: " << this->m_remaining_monoms);
}

template <std::size_t DIM>
bool BRY::UniformRandomFilter<DIM>::remove(const bry_int_t* exponent_vec) const {
    std::array<bry_int_t, DIM> idx;
    std::copy(exponent_vec, exponent_vec + DIM, idx.begin());
    return this->m_flags(idx);
}
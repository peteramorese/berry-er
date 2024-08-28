#pragma once

#include "MonomialFilter.h"

#include "berry/MultiIndex.h"
#include "berry/Operations.h"

#include <algorithm>

template <std::size_t DIM>
BRY::MonomialFilter<DIM>::MonomialFilter(bry_int_t barrier_deg) 
    : m_barrier_deg(barrier_deg)
    , m_flags(pow(barrier_deg + 1, DIM), false)
    , m_idx_map(pow(barrier_deg + 1, DIM))
{}

template <std::size_t DIM>
void BRY::MonomialFilter<DIM>::init() {
    m_remaining_monoms = 0;
    for (auto midx = mIdxW(DIM, m_barrier_deg + 1); !midx.last(); ++midx) {
        bool remove_monom = remove(midx.begin());
        m_flags[midx.inc().wrappedIdx()] = remove_monom;
        m_idx_map[midx.inc().wrappedIdx()] = m_remaining_monoms;
        m_remaining_monoms += !remove_monom;
        DEBUG("Idx: " << midx.inc().wrappedIdx() << " remove: " << remove_monom);
    }
    DEBUG("n remaining monoms: " << m_remaining_monoms);
}

template <std::size_t DIM>
BRY::bry_int_t BRY::MonomialFilter<DIM>::nRemainingMonoms() const {
    return m_remaining_monoms;
}

template <std::size_t DIM>
const std::vector<bool>& BRY::MonomialFilter<DIM>::flags() const {
    return m_flags;
}

template <std::size_t DIM>
const BRY::bry_int_t BRY::MonomialFilter<DIM>::newWrappedIdx(bry_int_t old_wrapped_idx) const {
    return m_idx_map[old_wrapped_idx];
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
bool BRY::OddSumFilter<DIM>::remove(const bry_int_t* exponent_vec) const {
    return (std::accumulate(exponent_vec, exponent_vec + DIM, 0) + 1) % 2 == 0;
}
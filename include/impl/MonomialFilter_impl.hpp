#pragma once

#include "MonomialFilter.h"

#include "berry/MultiIndex.h"
#include "berry/Operations.h"

#include <algorithm>

template <std::size_t DIM>
BRY::DiagDegFilter<DIM>::DiagDegFilter(bry_int_t barrier_deg)
    : m_barrier_deg(barrier_deg)
{}

template <std::size_t DIM>
bool BRY::DiagDegFilter<DIM>::remove(const bry_int_t* exponent_vec) {
    return std::accumulate(exponent_vec, exponent_vec + DIM, 0) > m_barrier_deg;
}

template <std::size_t DIM>
bool BRY::OddSumFilter<DIM>::remove(const bry_int_t* exponent_vec) {
    return (std::accumulate(exponent_vec, exponent_vec + DIM, 0) + 1) % 2 == 0;
}
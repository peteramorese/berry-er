#pragma once

#include "HyperRectangle.h"
#include "MonomialFilter.h"

#include <list>

namespace BRY {

template <std::size_t DIM>
static std::list<HyperRectangle<DIM>> makeRectBoundary(const HyperRectangle<DIM>& workspace, bry_float_t boundary_width, bry_int_t degree_increase = 0);

template <std::size_t DIM>
static void removeFilterFromPVector(Vector& polynomial_vector, const MonomialFilter<DIM>& filter);

} // namespace BRY

#include "impl/Tools_impl.hpp"
#pragma once

#include "HyperRectangle.h"

#include <list>

namespace BRY {

template <std::size_t DIM>
std::list<HyperRectangle<DIM>> makeRectBoundary(const HyperRectangle<DIM>& workspace, bry_float_t boundary_width, bry_int_t degree_increase = 0);

} // namespace BRY

#include "impl/Tools_impl.hpp"
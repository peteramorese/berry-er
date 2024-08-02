#pragma once

#include "HyperRectangle.h"

#include <list>

namespace BRY {

template <std::size_t DIM>
std::list<HyperRectangle<DIM>> makeRectBoundary(const HyperRectangle<DIM>& workspace, bry_float_t boundary_width);

    
} // namespace BRY

#include "impl/Tools_impl.hpp"
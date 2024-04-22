#pragma once

#include "HyperRectangle.h"

#include "berry/Polynomial.h"

#include <Eigen/Core>

namespace BRY {

namespace Transform {

template <std::size_t DIM>
static Eigen::MatrixXd affine(const HyperRectangle<DIM>& set, bry_deg_t m);

}

}

#include "impl/Transform_impl.hpp"
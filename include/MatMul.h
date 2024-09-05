#pragma once

#ifdef BRYR_USE_CUBLAS

#include <Eigen/Dense>

#include "berry/Types.h"

namespace BRY {

    Matrix multiplyCUDA(const Matrix& A, const Matrix& B);

}

// Use cuBLAS algorithms for matrix multiplication
#define BRY_MATMUL(A, B) BRY::multiplyCUDA(A, B)

#else

// Use regular Eigen algorithms for matrix multiplication
#define BRY_MATMUL(A, B) A * B

#endif

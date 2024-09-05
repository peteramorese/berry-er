#include "MatMul.h"

#include "lemon/Logging.h"

#ifdef BRYR_USE_CUBLAS

#include <exception>

//BRY::Matrix BRY::multiplyCUDA(const Matrix& A, const Matrix& B) {DEBUG("BEEPOS:LK");}

#include <cublas_v2.h>

BRY::Matrix BRY::multiplyCUDA(const Matrix& A, const Matrix& B) {
    if (A.cols() != B.rows()) {
        throw std::invalid_argument("Matrix dimension mismatch");
    }

    bry_int_t m = A.rows();
    bry_int_t k = A.cols();
    bry_int_t n = B.cols();

    BRY::Matrix result(m, n);
    DEBUG("b4 malloc");
    double *d_1, *d_2, *d_result;
    cudaMalloc((void**)&d_1, m * k * sizeof(bry_float_t));
    cudaMalloc((void**)&d_2, k * n * sizeof(bry_float_t));
    cudaMalloc((void**)&d_result, m * n * sizeof(bry_float_t));
    DEBUG("af malloc");

    // Copy matrices from Eigen to GPU
    DEBUG("b4 copy");
    cudaMemcpy(d_1, A.data(), m * k * sizeof(bry_float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_2, B.data(), k * n * sizeof(bry_float_t), cudaMemcpyHostToDevice);
    DEBUG("af copy");

    cublasStatus_t status;

    // Initialize cuBLAS
    cublasHandle_t handle;
    status = cublasCreate(&handle);
    if (status != CUBLAS_STATUS_SUCCESS) {
        ERROR("cuBLAS create handle failed with status " << status);
        cudaFree(d_1);
        cudaFree(d_2);
        cudaFree(d_result);
        throw std::runtime_error("cuBLAS error");
    }

    // Matrix multiplication: C = alpha * A * B + beta * C
    DEBUG("b4 gemm");
    bry_float_t alpha = 1.0;
    bry_float_t beta = 0.0;
    status = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, d_1, m, d_2, k, &beta, d_result, m);
    if (status != CUBLAS_STATUS_SUCCESS) {
        ERROR("cuBLAS gemm failed with status " << status);
        cudaFree(d_1);
        cudaFree(d_2);
        cudaFree(d_result);
        throw std::runtime_error("cuBLAS error");
    }
    DEBUG("af gemm");

    // Copy result from GPU to Eigen
    DEBUG("b4 copy out");
    cudaMemcpy(result.data(), d_result, m * n * sizeof(bry_float_t), cudaMemcpyDeviceToHost);
    DEBUG("af copy out");

    // Cleanup
    cublasDestroy(handle);
    cudaFree(d_1);
    cudaFree(d_2);
    cudaFree(d_result);

    return result;
}

#endif


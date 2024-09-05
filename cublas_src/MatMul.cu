#include "MatMul.h"

#include "lemon/Logging.h"

#ifdef BRYR_USE_CUBLAS

#include <exception>

#include <cublas_v2.h>

BRY::Matrix BRY::multiplyCUDA(const Matrix& A, const Matrix& B) {
    if (A.cols() != B.rows()) {
        throw std::invalid_argument("Matrix dimension mismatch");
    }

    bry_int_t m = A.rows();
    bry_int_t k = A.cols();
    bry_int_t n = B.cols();

    BRY::Matrix result(m, n);
    double *d_1, *d_2, *d_result;
    cudaMalloc((void**)&d_1, m * k * sizeof(bry_float_t));
    cudaMalloc((void**)&d_2, k * n * sizeof(bry_float_t));
    cudaMalloc((void**)&d_result, m * n * sizeof(bry_float_t));

    // Copy matrices to GPU
    cudaMemcpy(d_1, A.data(), m * k * sizeof(bry_float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_2, B.data(), k * n * sizeof(bry_float_t), cudaMemcpyHostToDevice);

    cublasStatus_t status;

    // Initialize cuBLAS
    cublasHandle_t handle;
    status = cublasCreate(&handle);
    if (status != CUBLAS_STATUS_SUCCESS) {
        ERROR("cuBLAS create handle failed with status " << status);
        cudaFree(d_1);
        cudaFree(d_2);
        cudaFree(d_result);
        throw std::runtime_error("cuBLAS initialization error");
    }

    // Matrix multiplication: C = alpha * A * B + beta * C
    bry_float_t alpha = 1.0;
    bry_float_t beta = 0.0;
    status = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, d_1, m, d_2, k, &beta, d_result, m);
    if (status != CUBLAS_STATUS_SUCCESS) {
        ERROR("cuBLAS gemm failed with status " << status);
        cudaFree(d_1);
        cudaFree(d_2);
        cudaFree(d_result);
        throw std::runtime_error("cuBLAS gemm error");
    }

    // Copy result from GPU to Eigen
    cudaMemcpy(result.data(), d_result, m * n * sizeof(bry_float_t), cudaMemcpyDeviceToHost);

    // Cleanup
    cublasDestroy(handle);
    cudaFree(d_1);
    cudaFree(d_2);
    cudaFree(d_result);

    return result;
}

#endif


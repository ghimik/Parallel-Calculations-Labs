#include "../include/matrix.hpp"
#include "../include/timer.hpp"
#include <cuda_runtime.h>
#include <iostream>
#include <cmath>

#define CUDA_CHECK(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line)
{
    if (code != cudaSuccess) 
    {
        std::cerr << "CUDA Error: " << cudaGetErrorString(code) << " " << file << ":" << line << std::endl;
        exit(code);
    }
}

__global__ void normalize_pivot_row(double* aug, double* res, int n, int k, double pivot) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j < n) {
        aug[k*n + j] /= pivot;
        res[k*n + j] /= pivot;
    }
}

__global__ void eliminate_rows(double* aug, double* res, int n, int k) {
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < n && j < n && i != k) {
        double factor = aug[i*n + k];
        aug[i*n + j] -= factor * aug[k*n + j];
        res[i*n + j] -= factor * res[k*n + j];
    }
}

Matrix inverse_cuda(const Matrix& A) {
    int n = A.n;
    Matrix aug = A;
    Matrix result = Matrix::identity(n);

    double* d_aug;
    double* d_res;

    CUDA_CHECK(cudaMalloc(&d_aug, n*n*sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_res, n*n*sizeof(double)));

    CUDA_CHECK(cudaMemcpy(d_aug, aug.data.data(), n*n*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_res, result.data.data(), n*n*sizeof(double), cudaMemcpyHostToDevice));

    for (int k = 0; k < n; k++) {
        int maxRow = k;
        double maxVal = std::abs(aug(k,k));
        for (int i = k+1; i < n; i++) {
            if (std::abs(aug(i,k)) > maxVal) {
                maxVal = std::abs(aug(i,k));
                maxRow = i;
            }
        }
        if (maxVal < 1e-12) throw std::runtime_error("Matrix is singular");

        if (maxRow != k) {
            for (int j = 0; j < n; j++) {
                std::swap(aug(k,j), aug(maxRow,j));
                std::swap(result(k,j), result(maxRow,j));
            }
            CUDA_CHECK(cudaMemcpy(d_aug, aug.data.data(), n*n*sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_res, result.data.data(), n*n*sizeof(double), cudaMemcpyHostToDevice));
        }

        double pivot = aug(k,k);

        int blockSize = 256;
        int gridSize = (n + blockSize - 1)/blockSize;
        normalize_pivot_row<<<gridSize, blockSize>>>(d_aug, d_res, n, k, pivot);
        CUDA_CHECK(cudaDeviceSynchronize());

        dim3 block(16,16);
        dim3 grid((n+block.x-1)/block.x, (n+block.y-1)/block.y);
        eliminate_rows<<<grid, block>>>(d_aug, d_res, n, k);
        CUDA_CHECK(cudaDeviceSynchronize());

        CUDA_CHECK(cudaMemcpy(&aug(k,0), &d_aug[k*n], n*sizeof(double), cudaMemcpyDeviceToHost));
    }

    CUDA_CHECK(cudaMemcpy(result.data.data(), d_res, n*n*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaFree(d_aug));
    CUDA_CHECK(cudaFree(d_res));

    return result;
}

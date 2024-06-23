// #include "chap4.hpp"

#include "gpu.hpp"
#include <cstdlib>

#include <cuda_runtime.h>

// c4v0 //////////////////////////////////////////////////////////////

__global__ void kernel_c4v0(int n, const float *data, float *result)
{
    // 1. Convert block/thread coordinate to global coordinate
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= n || j >= n)
        return;

    // 2. Start calculating values
    float v = HUGE_VALF;
    for (int k = 0; k < n; ++k)
    {
        v = min(v, data[n * i + k] + data[n * k + j]);
    }

    result[n * i + j] = v;
}

void GPU::c4v0(int n, const float *data, float *result)
{
    printf("Version: c4v0-own\n");

    int nBytes = n * n * sizeof(float);

    // 1. Allocate memory on GPU & copy date CPU -> GPU
    float *dGPU = nullptr, *rGPU = nullptr;
    cudaMalloc((void **)&dGPU, nBytes);
    cudaMalloc((void **)&rGPU, nBytes);

    cudaMemcpy(dGPU, data, nBytes, cudaMemcpyHostToDevice);

    // 2. Define block size & run kernel
    int nBlocks = (n + 16 - 1) / 16;
    dim3 dimThreads(16, 16);
    dim3 dimBlocks(nBlocks, nBlocks);

    kernel_c4v0<<<dimBlocks, dimThreads>>>(n, dGPU, rGPU);

    // 3. Copy data from GPU -> CPU
    cudaMemcpy(result, rGPU, nBytes, cudaMemcpyDeviceToHost);

    // 4. Release mem in GPU
    cudaFree(dGPU);
    cudaFree(rGPU);
}

// c4v1 //////////////////////////////////////////////////////////////

__global__ void kernel_c4v1(int n, const float *data, float *result)
{
    // 1. Convert block/thread coordinate to global coordinate
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= n || j >= n)
        return;

    // 2. Start calculating values
    float v = HUGE_VALF;
    for (int k = 0; k < n; ++k)
    {
        v = min(v, data[n * j + k] + data[n * k + i]);
    }

    result[n * i + j] = v;
}

void GPU::c4v1(int n, const float *data, float *result)
{
    printf("Version: c4v1-own\n");

    int nBytes = n * n * sizeof(float);

    // 1. Allocate memory on GPU & copy date CPU -> GPU
    float *dGPU = nullptr, *rGPU = nullptr;
    cudaMalloc((void **)&dGPU, nBytes);
    cudaMalloc((void **)&rGPU, nBytes);

    cudaMemcpy(dGPU, data, nBytes, cudaMemcpyHostToDevice);

    // 2. Define block size & run kernel
    int nBlocks = (n + 16 - 1) / 16;
    dim3 dimThreads(16, 16);
    dim3 dimBlocks(nBlocks, nBlocks);

    kernel_c4v0<<<dimBlocks, dimThreads>>>(n, dGPU, rGPU);

    // 3. Copy data from GPU -> CPU
    cudaMemcpy(result, rGPU, nBytes, cudaMemcpyDeviceToHost);

    // 4. Release mem in GPU
    cudaFree(dGPU);
    cudaFree(rGPU);
}

// c4v2 //////////////////////////////////////

void mykernel(float *r, const float *d, int n, int nn)
{
    int tX = threadIdx.x;
    int tY = threadIdx.y;
    int bX = blockIdx.x;
    int bY = blockIdx.y;

    const float *t = d + nn * nn;

    float v[8][8];
    for (int i_ = 0; i_ < 8; ++i_)
    {
        for (int j_ = 0; j_ < 8; ++j_)
        {
            v[i_][j_] = HUGE_VALF;
        }
    }
    for (int k = 0; k < n; ++k)
    {
        float x[8];
        float y[8];
        for (int i_ = 0; i_ < 8; ++i_)
        {
            int i = bX * 64 + i_ * 8 + tX;
            x[i_] = t[nn * k + i];
        }
        for (int j_ = 0; j_ < 8; ++j_)
        {
            int j = bY * 64 + j_ * 8 + tY;
            y[j_] = d[nn * k + j];
        }
        for (int i_ = 0; i_ < 8; ++i_)
        {
            for (int j_ = 0; j_ < 8; ++j_)
            {
                v[i_][j_] = min(v[i_][j_], x[i_] + y[j_]);
            }
        }
    }
    for (int i_ = 0; i_ < 8; ++i_)
    {
        for (int j_ = 0; j_ < 8; ++j_)
        {
            int i = bX * 64 + i_ * 8 + tX;
            int j = bY * 64 + j_ * 8 + tY;
            if (i < n && j < n)
            {
                r[n * i + j] = v[i_][j_];
            }
        }
    }
}
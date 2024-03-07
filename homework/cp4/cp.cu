#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

static inline void check(cudaError_t err, const char *context)
{
    if (err != cudaSuccess)
    {
        std::cerr << "CUDA error: " << context << ": "
                  << cudaGetErrorString(err) << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

#define CHECK(x) check(x, #x)

static inline int divup(int a, int b)
{
    return (a + b - 1) / b;
}

__global__ void kernel(float *d, float *r, int ny, int nx)
{
    // Convert to global coordinate
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (j < i || i >= ny || j >= ny)
    {
        return;
    }
    else if (j == i)
    {
        r[ny * i + j] = 1;
    }
    else
    {
        float cor = 0;
        for (int k = 0; k < nx; ++k)
            cor += d[i * nx + k] * d[j * nx + k];

        r[i * ny + j] = cor;
    }
}

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/
void correlate(int ny, int nx, const float *data, float *result)
{
    float *norm = new float[ny * nx];

    // 1. row-wise 0-mean normalization

    // #pragma omp parallel for
    for (int i = 0; i < ny; ++i)
    {
        float mean = 0;
        for (int j = 0; j < nx; ++j)
            // #pragma omp critical
            mean += data[i * nx + j];

        mean /= nx;

        for (int j = 0; j < nx; ++j)
            norm[i * nx + j] = data[i * nx + j] - mean;
    }

    // 2. row-wise square-sum normalization

    // #pragma omp parallel for
    for (int i = 0; i < ny; ++i)
    {
        float sq_sum = 0;
        for (int j = 0; j < nx; ++j)
            // #pragma omp critical
            sq_sum += pow(norm[i * nx + j], 2);
        sq_sum = sqrt(sq_sum);
        for (int j = 0; j < nx; ++j)
            norm[i * nx + j] = norm[i * nx + j] / sq_sum;
    }

    // 3. upper-triangular matmul

    // 3.1. Allocate device memory
    size_t size_norm = ny * nx * sizeof(float);
    size_t size_result = ny * ny * sizeof(float);

    float *d_norm = nullptr, *d_r = nullptr;

    CHECK(cudaMalloc(&d_norm, size_norm));
    CHECK(cudaMalloc(&d_r, size_result));

    // 3.2. Move data from host -> device
    const int BLOCK_SIZE = 16;
    int grid_size = divup(ny, BLOCK_SIZE);

    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    dim3 dimGrid(grid_size, grid_size);

    cudaMemcpy(d_norm, norm, size_norm, cudaMemcpyHostToDevice);
    cudaMemcpy(d_r, result, size_result, cudaMemcpyHostToDevice);

    // 3.3. Execute kernel
    kernel<<<dimGrid, dimBlock>>>(d_norm, d_r, ny, nx);

    cudaDeviceSynchronize();

    // 3.4. Move data from device back to host
    cudaMemcpy(result, d_r, size_result, cudaMemcpyDeviceToHost);

    CHECK(cudaFree(d_r));
    CHECK(cudaFree(d_norm));
}

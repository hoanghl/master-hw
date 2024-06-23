#include <cuda.h>
#include <cuda_runtime_api.h>
#include <iostream>

using namespace std;

#ifndef GPU_HPP
#define GPU_HPP

// static inline void check(cudaError_t err, const char *context)
// {
//     if (err != cudaSuccess)
//     {
//         std::cerr << "CUDA error: " << context << ": "
//                   << cudaGetErrorString(err) << std::endl;
//         std::exit(EXIT_FAILURE);
//     }
// }

static inline int divup(int a, int b)
{
    return (a + b - 1) / b;
}

// #define CHECK(x) check(x, #x)

namespace GPU
{
    void c4v0(int n, const float *data, float *result);
    void c4v1(int n, const float *data, float *result);
    // void c4v2(int n, const float *data, float *result);
    // void c4v3(int n, const float *data, float *result);
}

#endif

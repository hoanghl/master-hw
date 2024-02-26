#include <cmath>
#include <vector>
#include <x86intrin.h>
using namespace std;

constexpr int N_PER_PACK = 4;

typedef double double4 __attribute__((vector_size(N_PER_PACK * sizeof(double))));

const double4 vZeros = {
    0.,
    0.,
    0.,
    0.,
};

// static inline double4 swap4(double4 x) { return _mm256_permute2f128_ps(x, x, 0b00000001); }
static inline double4 swap2(double4 x) { return _mm256_permute4x64_pd(x, 0b01001110); }
static inline double4 swap1(double4 x) { return _mm256_permute4x64_pd(x, 0b10110001); }

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
    vector<double> norm(ny * nx, 0.);

    // 1. row-wise 0-mean normalization

#pragma omp parallel for
    for (int i = 0; i < ny; ++i)
    {
        double mean = 0;
        for (int j = 0; j < nx; ++j)
            mean += data[i * nx + j];

        mean /= nx;

        for (int j = 0; j < nx; ++j)
            norm[i * nx + j] = data[i * nx + j] - mean;
    }

    // 2. row-wise square-sum normalization

#pragma omp parallel for
    for (int i = 0; i < ny; ++i)
    {
        double sq_sum = 0;
        for (int j = 0; j < nx; ++j)
        {
            sq_sum += norm[i * nx + j] * norm[i * nx + j];
        }
        sq_sum = sqrt(sq_sum);
        for (int j = 0; j < nx; ++j)
            norm[i * nx + j] = norm[i * nx + j] / sq_sum;
    }

    // 3. upper-triangular matmul
    // 3.1. Expand vector 'norm'
    constexpr int bs = 4;
    int nPacks = (ny + N_PER_PACK - 1) / N_PER_PACK;

    // 3.2. Move 'norm' to new padded
    vector<double4> vnorm(nPacks * nx, vZeros);
#pragma omp parallel for
    for (int i = 0; i < ny; ++i)
        for (int j = 0; j < nx; ++j)
            vnorm[nx * (i / N_PER_PACK) + j][i % N_PER_PACK] = norm[nx * i + j];

// 3.3. Start calculating
#pragma omp parallel for
    for (int i = 0; i < ny; ++i)
        result[ny * i + i] = 1.0;

#pragma omp parallel for
    for (int i1 = 0; i1 < nPacks; ++i1)
        for (int i2 = i1; i2 < nPacks; ++i2)
        {
            double4 v0 = vZeros;
            double4 v1 = vZeros;
            double4 v2 = vZeros;
            double4 v3 = vZeros;
            // double4 v4 = vZeros;
            // double4 v5 = vZeros;
            // double4 v6 = vZeros;
            // double4 v7 = vZeros;

            // Looping
            for (int k = 0; k < nx; ++k)
            {
                double4 a0 = vnorm[nx * i1 + k], b0 = vnorm[nx * i2 + k];
                double4 a2 = swap2(a0);
                // double4 a4 = swap4(a0);
                // double4 a6 = swap2(a4);
                double4 b1 = swap1(b0);

                v0 += a0 * b0;
                v1 += a0 * b1;
                v2 += a2 * b0;
                v3 += a2 * b1;
                // v4 += a4 + b0;
                // v5 += a4 + b1;
                // v6 += a6 + b0;
                // v7 += a6 + b1;
            }

            double4 v[] = {v0, v1, v2, v3};
            for (int j = 1; j < bs; j += 2)
                v[j] = swap1(v[j]);

            // Assign back to 'result'
            for (int i_bs = 0; i_bs < bs; ++i_bs)
                for (int j_bs = 0; j_bs < bs; ++j_bs)
                {
                    int i = i1 * N_PER_PACK + i_bs;
                    int j = i2 * N_PER_PACK + j_bs;

                    if (i < ny && j < ny && i < j)
                        result[ny * i + j] = (float)v[i_bs ^ j_bs][j_bs];
                }
        }
}
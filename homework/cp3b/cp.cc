#include <algorithm>
#include <cmath>
#include <tuple>
#include <vector>
#include <x86intrin.h>
using namespace std;

constexpr int N_PER_PACK = 16;

typedef float float16 __attribute__((vector_size(N_PER_PACK * sizeof(float))));

const float16 vZeros = {
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
};

static inline float16 swap1(float16 x) { return _mm512_permute_ps(x, 0b10110001); }
static inline float16 swap2(float16 x) { return _mm512_permute_ps(x, 0b01001110); }
static inline float16 swap4(float16 x) { return _mm512_shuffle_f32x4(x, x, 0b10110001); }
static inline float16 swap8(float16 x) { return _mm512_shuffle_f32x4(x, x, 0b01001110); }

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
    vector<float> norm(ny * nx, 0.);

    // 1. row-wise 0-mean normalization

#pragma omp parallel for
    for (int i = 0; i < ny; ++i)
    {
        float mean = 0;
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
        float sq_sum = 0;
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
    constexpr int bs = N_PER_PACK;
    int nPacks = (ny + N_PER_PACK - 1) / N_PER_PACK;

    // 3.2. Move 'norm' to new padded
    vector<float16> vnorm(nPacks * nx, vZeros);
#pragma omp parallel for
    for (int i = 0; i < ny; ++i)
        for (int j = 0; j < nx; ++j)
            vnorm[nx * (i / N_PER_PACK) + j][i % N_PER_PACK] = norm[nx * i + j];

// 3.3. Start calculating
#pragma omp parallel for
    for (int i = 0; i < ny; ++i)
        result[ny * i + i] = 1.0;

        //     std::vector<std::tuple<int, int, int>> rows(nPacks * nPacks);
        // #pragma omp parallel for
        //     for (int i1 = 0; i1 < nPacks; ++i1)
        //         for (int i2 = 0; i2 < nPacks; ++i2)
        //         {
        //             int ija = _pdep_u32(i1, 0x55555555) | _pdep_u32(i2, 0xAAAAAAAA);
        //             rows[i1 * nPacks + i2] = std::make_tuple(ija, i1, i2);
        //         }
        // std::sort(rows.begin(), rows.end());

#pragma omp parallel for
    for (int i1 = 0; i1 < nPacks; ++i1)
        for (int i2 = i1; i2 < nPacks; ++i2)
        {
            // for (int ij = 0; ij < nPacks * nPacks; ++ij)
            // {
            //     int i1 = std::get<1>(rows[ij]);
            //     int i2 = std::get<2>(rows[ij]);

            float16 v0 = vZeros;
            float16 v1 = vZeros;
            float16 v2 = vZeros;
            float16 v3 = vZeros;
            float16 v4 = vZeros;
            float16 v5 = vZeros;
            float16 v6 = vZeros;
            float16 v7 = vZeros;
            float16 v8 = vZeros;
            float16 v9 = vZeros;
            float16 v10 = vZeros;
            float16 v11 = vZeros;
            float16 v12 = vZeros;
            float16 v13 = vZeros;
            float16 v14 = vZeros;
            float16 v15 = vZeros;

            // Looping
            for (int k = 0; k < nx; ++k)
            {
                // constexpr int PF = 20;
                // __builtin_prefetch(&vnorm[nx * i2 + k + PF]);
                float16 a0 = vnorm[nx * i1 + k], b0 = vnorm[nx * i2 + k];
                float16 b1 = swap1(b0);
                float16 a2 = swap2(a0);
                float16 a4 = swap4(a0);
                float16 a6 = swap2(a4);
                float16 a8 = swap8(a0);
                float16 a10 = swap2(a8);
                float16 a12 = swap4(a8);
                float16 a14 = swap2(a12);

                v0 += a0 * b0;
                v1 += a0 * b1;
                v2 += a2 * b0;
                v3 += a2 * b1;
                v4 += a4 * b0;
                v5 += a4 * b1;
                v6 += a6 * b0;
                v7 += a6 * b1;
                v8 += a8 * b0;
                v9 += a8 * b1;
                v10 += a10 * b0;
                v11 += a10 * b1;
                v12 += a12 * b0;
                v13 += a12 * b1;
                v14 += a14 * b0;
                v15 += a14 * b1;
            }

            float16 v[] = {
                v0,
                v1,
                v2,
                v3,
                v4,
                v5,
                v6,
                v7,
                v8,
                v9,
                v10,
                v11,
                v12,
                v13,
                v14,
                v15,
            };
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
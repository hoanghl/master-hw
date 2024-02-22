// #include <avxintrin.h>
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>
#include <x86intrin.h>

#include "chap2.hpp"

using namespace std;
using namespace std::chrono;

typedef float float8_t __attribute__((vector_size(8 * sizeof(float))));

const float inf = numeric_limits<float>::infinity();
const float8_t vfloat8_inf = {
    inf,
    inf,
    inf,
    inf,
    inf,
    inf,
    inf,
    inf,
};

static inline float8_t swap4(float8_t x) { return _mm256_permute2f128_ps(x, x, 0b00000001); }
static inline float8_t swap2(float8_t x) { return _mm256_permute_ps(x, 0b01001110); }
static inline float8_t swap1(float8_t x) { return _mm256_permute_ps(x, 0b10110001); }

static inline float8_t minVec(float8_t x, float8_t y)
{
    return x < y ? x : y;
}
static inline float minInternal(float8_t v)
{
    float min = v[0];
    for (int i = 1; i < 8; ++i)
        if (min > v[i])
            min = v[i];

    return min;
}

void step(float *r, const float *d_, int n)
{
    constexpr int nPerPack = 8;
    constexpr int bsize = 8;

    // Calculate necessary things
    int nPacks = (n + nPerPack - 1) / nPerPack;

    // Create new arrays and move data to these new arrays
    vector<float8_t> d(nPacks * n, vfloat8_inf), dT(nPacks * n, vfloat8_inf);
#pragma omp parallel for
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        {
            d[(i / nPerPack) * n + j][i % nPerPack] = d_[n * i + j];
            dT[i * nPacks + j / nPerPack][j % nPerPack] = d_[n * i + j];
        }

// Start calculating
#pragma omp parallel for
    for (int i1 = 0; i1 < nPacks; ++i1)
        for (int i2 = 0; i2 < nPacks; i2++)
        {
            float8_t v0 = vfloat8_inf;
            float8_t v1 = vfloat8_inf;
            float8_t v2 = vfloat8_inf;
            float8_t v3 = vfloat8_inf;
            float8_t v4 = vfloat8_inf;
            float8_t v5 = vfloat8_inf;
            float8_t v6 = vfloat8_inf;
            float8_t v7 = vfloat8_inf;

            for (int k = 0; k < n; ++k)
            {
                float8_t a0 = d[i1 * nPacks + k];
                float8_t b0 = dT[i2 * nPacks + k];

                float8_t b1 = swap1(b0);
                float8_t a2 = swap2(a0);
                float8_t a4 = swap4(a0);
                float8_t a6 = swap4(a2);

                v0 = minVec(v0, a0 + b0);
                v1 = minVec(v1, a0 + b1);
                v2 = minVec(v2, a2 + b0);
                v3 = minVec(v3, a2 + b1);
                v4 = minVec(v4, a4 + b0);
                v5 = minVec(v5, a4 + b1);
                v6 = minVec(v6, a6 + b0);
                v7 = minVec(v7, a6 + b1);
            }

            // Do post-processing with vectors v
            float8_t vs[8] = {v0, v1, v2, v3, v4, v5, v6, v7};
            for (int i = 1; i < 8; i += 2)
                vs[i] = swap1(vs[i]);

            // Assign back to r
            for (int i = 0; i < 8; ++i)
                for (int j = 0; j < 8; ++j)
                {
                    int i_r = i1 * nPerPack + i;
                    int j_r = i2 * nPerPack + j;
                    if (i_r < n && j_r < n)
                        r[n * i_r + j_r] = vs[i ^ j][j];
                }
        }
}

void step_reference(float *r, const float *d_, int n)
{
    // vectors per input column
    int na = (n + 8 - 1) / 8;

    // input data, padded, converted to vectors
    std::vector<float8_t> vd(na * n);
    // input data, transposed, padded, converted to vectors
    std::vector<float8_t> vt(na * n);

#pragma omp parallel for
    for (int ja = 0; ja < na; ++ja)
    {
        for (int i = 0; i < n; ++i)
        {
            for (int jb = 0; jb < 8; ++jb)
            {
                int j = ja * 8 + jb;
                vd[n * ja + i][jb] = j < n ? d_[n * j + i] : inf;
                vt[n * ja + i][jb] = j < n ? d_[n * i + j] : inf;
            }
        }
    }

#pragma omp parallel for
    for (int ia = 0; ia < na; ++ia)
    {
        for (int ja = 0; ja < na; ++ja)
        {
            float8_t vv000 = vfloat8_inf;
            float8_t vv001 = vfloat8_inf;
            float8_t vv010 = vfloat8_inf;
            float8_t vv011 = vfloat8_inf;
            float8_t vv100 = vfloat8_inf;
            float8_t vv101 = vfloat8_inf;
            float8_t vv110 = vfloat8_inf;
            float8_t vv111 = vfloat8_inf;
            for (int k = 0; k < n; ++k)
            {
                float8_t a000 = vd[n * ia + k];
                float8_t b000 = vt[n * ja + k];
                float8_t a100 = swap4(a000);
                float8_t a010 = swap2(a000);
                float8_t a110 = swap2(a100);
                float8_t b001 = swap1(b000);
                vv000 = minVec(vv000, a000 + b000);
                vv001 = minVec(vv001, a000 + b001);
                vv010 = minVec(vv010, a010 + b000);
                vv011 = minVec(vv011, a010 + b001);
                vv100 = minVec(vv100, a100 + b000);
                vv101 = minVec(vv101, a100 + b001);
                vv110 = minVec(vv110, a110 + b000);
                vv111 = minVec(vv111, a110 + b001);
            }
            float8_t vv[8] = {vv000, vv001, vv010, vv011, vv100, vv101, vv110, vv111};
            for (int kb = 1; kb < 8; kb += 2)
            {
                vv[kb] = swap1(vv[kb]);
            }
            for (int jb = 0; jb < 8; ++jb)
            {
                for (int ib = 0; ib < 8; ++ib)
                {
                    int i = ib + ia * 8;
                    int j = jb + ja * 8;
                    if (j < n && i < n)
                    {
                        r[n * i + j] = vv[ib ^ jb][jb];
                    }
                }
            }
        }
    }
}

int main(int argc, char const *argv[])
{
    chap2_result *result = readFile(STD_FILENAME);
    cout << "1. Finish reading data" << endl;
    float *r = new float[result->n * result->n];

    cout << "2. Step" << endl;
    auto t1 = high_resolution_clock::now();
    step(r, result->d, result->n);
    auto t2 = high_resolution_clock::now();

    duration<double, std::milli> ms_double = t2 - t1;
    cout << "Running time: " << ms_double.count() << endl;

    return 0;
}

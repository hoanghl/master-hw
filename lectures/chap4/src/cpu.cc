#include "chap4.hpp"

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <math.h>
#include <tuple>
#include <vector>
#include <x86intrin.h>

using namespace std;

const float inf = numeric_limits<float>::infinity();
typedef float float8_t __attribute__((vector_size(8 * sizeof(float))));

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

void CPU::baseline(int n, const float *data, float *result)
{

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        {
            float v = numeric_limits<float>::infinity();
            for (int k = 0; k < n; ++k)
            {
                v = min(v, data[n * i + k] + data[n * k + j]);
            }

            result[n * i + j] = v;
        }
}

void CPU::c2v7(int n, const float *data, float *result)
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
                vd[n * ja + i][jb] = j < n ? data[n * j + i] : inf;
                vt[n * ja + i][jb] = j < n ? data[n * i + j] : inf;
            }
        }
    }

    std::vector<std::tuple<int, int, int>> rows(na * na);

    for (int ia = 0; ia < na; ++ia)
    {
        for (int ja = 0; ja < na; ++ja)
        {
            int ija = _pdep_u32(ia, 0x55555555) | _pdep_u32(ja, 0xAAAAAAAA);
            rows[ia * na + ja] = std::make_tuple(ija, ia, ja);
        }
    }
    std::sort(rows.begin(), rows.end());

#pragma omp parallel for schedule(dynamic, na)
    for (int i = 0; i < na * na; ++i)
    {
        int ia = std::get<1>(rows[i]);
        int ja = std::get<2>(rows[i]);

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
            constexpr int PF = 20;
            __builtin_prefetch(&vd[n * ia + k + PF]);
            __builtin_prefetch(&vt[n * ja + k + PF]);

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
                    result[n * i + j] = vv[ib ^ jb][jb];
                }
            }
        }
    }
}
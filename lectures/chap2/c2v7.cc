
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <tuple>
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

    std::vector<std::tuple<int, int, int>> rows(na * na);
#pragma omp parallel for
    for (int ia = 0; ia < na; ++ia)
    {
        for (int ja = 0; ja < na; ++ja)
        {
            int ija = _pdep_u32(ia, 0x55555555) | _pdep_u32(ja, 0xAAAAAAAA);
            rows[ia * na + ja] = std::make_tuple(ija, ia, ja);
        }
    }
    std::sort(rows.begin(), rows.end());

    vector<float> rTmp(n * n, 0.);
    // for (int kStart = 0; kStart < n; kStart += 8)
    // {

#pragma omp parallel for
    for (int i = 0; i < na * na; ++i)
    {
        int ia = std::get<1>(rows[i]);
        int ja = std::get<2>(rows[i]);
        // for (int kStart = 0; kStart < n; kStart += 32)
        //     for (int ia = 0; ia < na; ++ia)
        //     {
        //         for (int ja = 0; ja < na; ++ja)
        //         {

        float8_t vv000 = vfloat8_inf;
        float8_t vv001 = vfloat8_inf;
        float8_t vv010 = vfloat8_inf;
        float8_t vv011 = vfloat8_inf;
        float8_t vv100 = vfloat8_inf;
        float8_t vv101 = vfloat8_inf;
        float8_t vv110 = vfloat8_inf;
        float8_t vv111 = vfloat8_inf;

        // int kEnd = std::min(n, kStart + 32);
        // for (int k = kStart; k < kEnd; ++k)
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
                    r[n * i + j] = vv[ib ^ jb][jb];
            }
        }
    }
}

// for (int i = 0; i < n * n; ++i)
//     r[i] = rTmp[i];
// }

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

#pragma omp parallel for
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
                    r[n * i + j] = vv[ib ^ jb][jb];
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

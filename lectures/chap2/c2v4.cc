#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

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
    const int nPerPack = 8;
    const int bsize = 3;

    // Calculate necessary things
    int nPacks = (n + nPerPack - 1) / nPerPack;
    int nRowsPadded = n;
    if (nRowsPadded % bsize != 0)
        nRowsPadded = ceil(n * 1.0 / bsize) * bsize;

    // Create new arrays and move data to these new arrays
    vector<float8_t> d(nRowsPadded * nPacks, vfloat8_inf), dT(nRowsPadded * nPacks, vfloat8_inf);
#pragma omp parallel for
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        {
            d[i * nPacks + j / nPerPack][j % nPerPack] = d_[n * i + j];
            dT[i * nPacks + j / nPerPack][j % nPerPack] = d_[n * j + i];
        }

// Start calculating
#pragma omp parallel for
    for (int i_d = 0; i_d < nRowsPadded; i_d += bsize)
        for (int i_dT = 0; i_dT < nRowsPadded; i_dT += bsize)
        {
            // Create v and initialize
            float8_t v[bsize][bsize];
            for (int i = 0; i < bsize; ++i)
                for (int j = 0; j < bsize; ++j)
                    v[i][j] = vfloat8_inf;

            // Start scanning each pack
            for (int k = 0; k < nPacks; ++k)
            {
                float8_t vd_1 = d[nPacks * (i_d + 0) + k];
                float8_t vd_2 = d[nPacks * (i_d + 1) + k];
                float8_t vd_3 = d[nPacks * (i_d + 2) + k];
                float8_t vdT_1 = dT[nPacks * (i_dT + 0) + k];
                float8_t vdT_2 = dT[nPacks * (i_dT + 1) + k];
                float8_t vdT_3 = dT[nPacks * (i_dT + 2) + k];
                v[0][0] = minVec(v[0][0], vd_1 + vdT_1);
                v[0][1] = minVec(v[0][1], vd_1 + vdT_2);
                v[0][2] = minVec(v[0][1], vd_1 + vdT_3);
                v[1][0] = minVec(v[1][0], vd_2 + vdT_1);
                v[1][1] = minVec(v[1][1], vd_2 + vdT_2);
                v[1][2] = minVec(v[1][2], vd_2 + vdT_3);
                v[2][0] = minVec(v[2][0], vd_3 + vdT_1);
                v[2][1] = minVec(v[2][1], vd_3 + vdT_2);
                v[2][2] = minVec(v[2][2], vd_3 + vdT_3);
            }

            // Collapse
            for (int i = 0; i < bsize; ++i)
                for (int j = 0; j < bsize; ++j)
                    if (i_d + i < n && i_dT + j < n)
                        r[n * (i_d + i) + i_dT + j] = minInternal(v[i][j]);
        }
}

void step_reference(float *r, const float *d_, int n)
{
    // elements per vector
    constexpr int nb = 8;
    // vectors per input row
    int na = (n + nb - 1) / nb;

    // block size
    constexpr int nd = 3;
    // how many blocks of rows
    int nc = (n + nd - 1) / nd;
    // number of rows after padding
    int ncd = nc * nd;

    // input data, padded, converted to vectors
    std::vector<float8_t> vd(ncd * na);
    // input data, transposed, padded, converted to vectors
    std::vector<float8_t> vt(ncd * na);

#pragma omp parallel for
    for (int j = 0; j < n; ++j)
    {
        for (int ka = 0; ka < na; ++ka)
        {
            for (int kb = 0; kb < nb; ++kb)
            {
                int i = ka * nb + kb;
                vd[na * j + ka][kb] = i < n ? d_[n * j + i] : inf;
                vt[na * j + ka][kb] = i < n ? d_[n * i + j] : inf;
            }
        }
    }
    for (int j = n; j < ncd; ++j)
    {
        for (int ka = 0; ka < na; ++ka)
        {
            for (int kb = 0; kb < nb; ++kb)
            {
                vd[na * j + ka][kb] = inf;
                vt[na * j + ka][kb] = inf;
            }
        }
    }

#pragma omp parallel for
    for (int ic = 0; ic < nc; ++ic)
    {
        for (int jc = 0; jc < nc; ++jc)
        {
            float8_t vv[nd][nd];
            for (int id = 0; id < nd; ++id)
            {
                for (int jd = 0; jd < nd; ++jd)
                {
                    vv[id][jd] = vfloat8_inf;
                }
            }
            for (int ka = 0; ka < na; ++ka)
            {
                float8_t y0 = vt[na * (jc * nd + 0) + ka];
                float8_t y1 = vt[na * (jc * nd + 1) + ka];
                float8_t y2 = vt[na * (jc * nd + 2) + ka];
                float8_t x0 = vd[na * (ic * nd + 0) + ka];
                float8_t x1 = vd[na * (ic * nd + 1) + ka];
                float8_t x2 = vd[na * (ic * nd + 2) + ka];
                vv[0][0] = minVec(vv[0][0], x0 + y0);
                vv[0][1] = minVec(vv[0][1], x0 + y1);
                vv[0][2] = minVec(vv[0][2], x0 + y2);
                vv[1][0] = minVec(vv[1][0], x1 + y0);
                vv[1][1] = minVec(vv[1][1], x1 + y1);
                vv[1][2] = minVec(vv[1][2], x1 + y2);
                vv[2][0] = minVec(vv[2][0], x2 + y0);
                vv[2][1] = minVec(vv[2][1], x2 + y1);
                vv[2][2] = minVec(vv[2][2], x2 + y2);
            }
            for (int id = 0; id < nd; ++id)
            {
                for (int jd = 0; jd < nd; ++jd)
                {
                    int i = ic * nd + id;
                    int j = jc * nd + jd;
                    if (i < n && j < n)
                    {
                        r[n * i + j] = minInternal(vv[id][jd]);
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

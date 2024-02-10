#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

#include "chap2.hpp"

using namespace std;
using namespace std::chrono;

typedef float float8_t __attribute__((vector_size(8 * sizeof(float))));

const float infty = numeric_limits<float>::infinity();
const float8_t f8_inf{
    infty,
    infty,
    infty,
    infty,
    infty,
    infty,
    infty,
    infty,
};

static inline float8_t minVector(float8_t x, float8_t y)
{
    return x < y ? x : y;
}

static inline float minInternal(float8_t x)
{
    float minVal = infty;
    for (int i = 0; i < 8; ++i)
        minVal = x[i] < minVal ? x[i] : minVal;

    return minVal;
}

void step(float *r, const float *d, int n)
{
    const int nPerPack = 8;
    int nPadded = n;
    if (nPadded % nPerPack != 0)
        nPadded = static_cast<int>(ceil(n * 1.0 / nPerPack)) * nPerPack;
    int nPacks = nPadded / nPerPack;

    // Move original data to pack
    vector<float8_t> v(n * nPacks, f8_inf), vT(n * nPacks, f8_inf);
#pragma omp parallel for
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        {
            v[nPacks * i + j / nPerPack][j % nPerPack] = d[n * i + j];
            vT[nPacks * i + j / nPerPack][j % nPerPack] = d[n * j + i];
        }

        // Start calculating
#pragma omp parallel for
    for (int i1 = 0; i1 < n; ++i1)
        for (int i2 = 0; i2 < n; ++i2)
        {
            float8_t vMin = f8_inf;
            for (int j = 0; j < nPacks; ++j)
            {
                float8_t x = v[nPacks * i1 + i2];
                float8_t y = vT[nPacks * i1 + i2];
                vMin = minVector(vMin, x + y);
            }

            r[n * i1 + i2] = minInternal(vMin);
        }
}

int main()
{
    chap2_result *result = readFile(STD_FILENAME);

    cout << "Finish reading data!" << endl;

    float *r = new float[result->n * result->n];

    cout << "Start stepping!" << endl;

    auto t1 = high_resolution_clock::now();
    step(r, result->d, result->n);
    auto t2 = high_resolution_clock::now();

    duration<double, std::milli> ms_double = t2 - t1;
    cout << "Running time: " << ms_double.count() << endl;

    return 0;
}
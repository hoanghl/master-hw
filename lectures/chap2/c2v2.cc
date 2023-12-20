#include <iostream>
#include <chrono>
#include <vector>

#include "chap2.hpp"

using namespace std;
using namespace std::chrono;

constexpr float infty = std::numeric_limits<float>::infinity();

void step(float *r, const float *d_, int n)
{
    constexpr int nb = 4;
    int na = (n + nb - 1) / nb;
    int nab = na * nb;

    // input data, padded
    std::vector<float> d(n * nab, infty);
    // input data, transposed, padded
    std::vector<float> t(n * nab, infty);

#pragma omp parallel for
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            d[nab * j + i] = d_[n * j + i];
            t[nab * j + i] = d_[n * i + j];
        }
    }

    auto t1 = high_resolution_clock::now();
#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            // vv[0] = result for k = 0, 4, 8, ...
            // vv[1] = result for k = 1, 5, 9, ...
            // vv[2] = result for k = 2, 6, 10, ...
            // vv[3] = result for k = 3, 7, 11, ...
            float vv[nb];
            for (int kb = 0; kb < nb; ++kb)
            {
                vv[kb] = infty;
            }
            for (int ka = 0; ka < na; ++ka)
            {
                for (int kb = 0; kb < nb; ++kb)
                {
                    float x = d[nab * i + ka * nb + kb];
                    float y = t[nab * j + ka * nb + kb];
                    float z = x + y;
                    vv[kb] = std::min(vv[kb], z);
                }
            }
            // v = result for k = 0, 1, 2, ...
            float v = infty;
            for (int kb = 0; kb < nb; ++kb)
            {
                v = std::min(vv[kb], v);
            }
            r[n * i + j] = v;
        }
    }
    auto t2 = high_resolution_clock::now();

    duration<double, std::milli> ms_double = t2 - t1;
    cout << "Running time: " << ms_double.count() << endl;
}

int main(int argc, char const *argv[])
{
    chap2_result *result = readFile(STD_FILENAME);
    cout << "Finish reading data!" << endl;

    float *r = new float[result->n * result->n];

    cout << "Start stepping!" << endl;

    step(r, result->d, result->n);

    return 0;
}

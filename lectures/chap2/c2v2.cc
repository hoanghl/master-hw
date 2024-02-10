#include "chap2.hpp"
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;
using namespace std::chrono;

void step(float *r, const float *d, int n)
{
    float inf = numeric_limits<float>::infinity();

    // 1. Calculate padding
    const uint32_t nCon = 4;
    int n_padded = ceil(n * 1.0 / nCon) * nCon;

    // cout << "n = " << n << " - "
    //      << "n_padded = " << n_padded << endl;

    // 2. Move to padded array and do transpose
    vector<float> d_padded(n * n_padded), d_t(n * n_padded);

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        {
            // if (n_padded * i + j > n * n_padded)
            // {
            //     cout << "Exceed: " << n_padded * i + j << endl;
            //     exit(1);
            // }
            d_t[n_padded * i + j] = d[n * j + i];
            d_padded[n_padded * i + j] = d[n * i + j];
        }

    // 3. Calculate
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            float v_arr[nCon];
            for (int v = 0; v < nCon; ++v)
                v_arr[v] = inf;

            for (int k = 0; k < ceil(n * 1.0 / nCon); ++k)
                for (int l = 0; l < nCon; ++l)
                {
                    float x = d[n_padded * i + nCon * k + l];
                    float y = d_t[n_padded * j + nCon * k + l];
                    v_arr[l] = min(v_arr[l], x + y);
                }

            r[n * i + j] = inf;
            for (int l = 0; l < nCon; ++l)
                r[n * i + j] = min(r[n * i + j], v_arr[l]);
        }
    }
}

int main(int argc, char const *argv[])
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

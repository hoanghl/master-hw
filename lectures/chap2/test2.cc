#include <iostream>
#include <chrono>
#include <vector>

#include "chap2.hpp"

using namespace std;
using namespace std::chrono;

constexpr float inf = std::numeric_limits<float>::infinity();

int main(int argc, char const *argv[])
{

    chap2_result *data = readFile(STD_FILENAME);

    constexpr int n_split = 4;
    int npad = 0;
    if (data->n % n_split != 0)
    {
        npad = ((data->n / n_split) + 1) * n_split - data->n;
    }
    int n_final = data->n + npad;
    vector<float> d(data->n * n_final, inf);
    vector<float> t(data->n * n_final, inf);

    vector<float> v_split(4, 0);
    int n_fold = n_final / n_split;

    // 1. Reallocate
    auto t1 = high_resolution_clock::now();
#pragma omp parallel for
    for (int i = 0; i < data->n; ++i)
        for (int j = 0; j < data->n; ++j)
        {
            d[n_final * i + j] = data->d[data->n * i + j];
            t[n_final * i + j] = data->d[data->n * j + i];
        }

        // 2. Start looping

#pragma omp parallel for
    for (int i = 0; i < data->n; ++i)
    {
        for (int j = 0; j < data->n; ++j)
        {
            for (int k1 = 0; k1 < n_split; ++k1)
            {
                for (int k2 = 0; k2 < n_fold; ++k2)
                {
                    float dist_ik = d[n_final * i + k2 * n_split + k1];
                    float dist_kj = t[n_final * i + k2 * n_split + k1];
                    v_split[k1] = std::min(v_split[k1], dist_ik + dist_kj);
                }
            }

            float v = inf;
            for (int k1 = 0; k1 < n_split; ++k1)
            {
                v = std::min(v_split[k1], v);
            }
        }
    }
    auto t2 = high_resolution_clock::now();

    duration<double, std::milli> ms_double = t2 - t1;
    cout << "Running time: " << ms_double.count() << endl;

    return 0;
}

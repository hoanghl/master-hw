#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <tuple>
#include <vector>

// #include <x86intrin.h>

// #include "exp.hpp"

using namespace std;
using namespace std::chrono;

// typedef double double4_t __attribute__((vector_size(4 * sizeof(double))));

// const double4_t v_double4_0 = {
//     0.,
//     0.,
//     0.,
//     0.,
// };

// static inline double sum_v(double4_t x) { return x[0] + x[1] + x[2]; }

static inline int getNumCells(int left, int top, int right, int bot)
{
    int nRowsInside = bot - top + 1;
    int nColsInside = right - left + 1;

    return nRowsInside * nColsInside;
}

int main(int argc, char const *argv[])
{
    int ny = 200, nx = 200;

    const int NYX = ny * nx;
    constexpr int N_VECS = 8;
    constexpr int N_PER_VECS = 4;

    int n = (NYX + N_VECS - 1) / N_VECS;
    int NYX_pad = n * N_VECS;
    const int NYX2 = NYX * NYX;

    auto t1 = high_resolution_clock::now();

    for (int size = 0; size < NYX; ++size)
        for (int t = 0; t < 1; ++t)
            for (int l = 0; l < 1; ++l)
                for (int height = 1; height <= min(ny - 1, size + t); ++height)
                {
                    // Determine b, r
                    int b, r;
                    if (size % height != 0)
                        continue;
                    else
                    {
                        int width = size / height;
                        r = l + width - 1;
                        b = t + height - 1;

                        if (r >= nx || b >= ny)
                            continue;
                    }
                }

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    printf("1. Running time: %8.4f\n", ms_double.count());
}

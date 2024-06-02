#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <tuple>
#include <vector>

#include <x86intrin.h>

#include "exp.hpp"

using namespace std;
using namespace std::chrono;

typedef double double4_t __attribute__((vector_size(4 * sizeof(double))));

const double4_t v_double4_0 = {
    0.,
    0.,
    0.,
    0.,
};

static inline double sum_v(double4_t x) { return x[0] + x[1] + x[2]; }

Result segment(int ny, int nx, const float *data)
{
    // NOTE: HoangLe [May-24]: In my code, bottom-r is in included mode whereas the requirement of the output is excluded mode

    const int NYX = ny * nx;
    constexpr int N_VECS = 8;
    constexpr int N_PER_VECS = 4;

    int n = (NYX + N_VECS - 1) / N_VECS;
    int NYX_pad = n * N_VECS;
    const int NYX2 = NYX * NYX;

    Result result{0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}};

    // 1. Pre-compute

    // 1.1. Define pre-computed vectors
    vector<double4_t> cumSum = vector<double4_t>(nx * ny, v_double4_0);
    // cumSum[0] = {data[0], data[1], data[2], 0.};
    cumSum[0][0] = data[0];
    cumSum[0][1] = data[1];
    cumSum[0][2] = data[2];

    // 1.2. Pre-compute
    auto t1 = high_resolution_clock::now();

    // NOTE: HoangLe [May-21]: Can use multithread, z-order here

    for (int br = 1; br < NYX; ++br) // Start at 2 because position [0, 0] is pre-assigned above
    {
        int r, b;
        Utils::xy2each(br, nx, r, b);
        int t = b - 1, l = r - 1;

        // Determine regions' location
        int locX = t < 0 || l < 0 ? -1 : l + nx * t;
        int locXY = t < 0 ? -1 : r + nx * t;
        int locXZ = l < 0 ? -1 : l + nx * b;
        int locCurrent = r + nx * b;

        double4_t X = locX == -1 ? v_double4_0 : cumSum[locX];
        double4_t XY = locXY == -1 ? v_double4_0 : cumSum[locXY];
        double4_t XZ = locXZ == -1 ? v_double4_0 : cumSum[locXZ];
        double4_t current = {data[0 + 3 * locCurrent], data[1 + 3 * locCurrent], data[2 + 3 * locCurrent], 0.};

        cumSum[locCurrent] = XY + XZ - X + current;
    }

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    printf("1. Running time: %8.4f\n", ms_double.count());

    // 2. Assign inner, outer, cost

    t1 = high_resolution_clock::now();

    vector<Result> bestResult = vector<Result>(NYX_pad);
    vector<double> costs = vector<double>(NYX, inf);

    // TODO: HoangLe [May-25]: Fill 3 above vectors with heavy optimizations

#pragma omp parallel for schedule(dynamic, ny)
    for (int size = 1; size < NYX; ++size)
    {
        double bestCostEach = inf;
        Result bestResultEach{0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}};
        for (int tlbr = 0; tlbr < NYX2; ++tlbr)
        {
            int tl, br, t, l, b, r;
            tl = tlbr / NYX;
            br = tlbr % NYX;
            t = tl / nx;
            l = tl % nx;
            b = br / nx;
            r = br % nx;

            if (b < t || r < l || t >= ny || b >= ny || l >= nx || r >= nx)
                continue;

            int nInside = Utils::getNumCells(l, t, r, b);
            int nOutside = Utils::getNumCellsOutside(l, t, r, b, nx, ny);

            if (nOutside == 0 || nInside != size)
                continue;

            double4_t v_nInside = _mm256_set1_pd(nInside);
            double4_t v_nOutside = _mm256_set1_pd(nOutside);

            double4_t innerEach = v_double4_0;
            double4_t outerEach = v_double4_0;

            // Calculate components for include-exclude principle
            int locX = l == 0 || t == 0 ? -1 : (l - 1) + nx * (t - 1);
            int locXY = t == 0 ? -1 : r + nx * (t - 1);
            int locXZ = l == 0 ? -1 : (l - 1) + nx * b;
            int locXYZW = r + nx * b;
            int locWhole = nx * ny - 1;

            double4_t X = locX == -1 ? v_double4_0 : cumSum[locX];
            double4_t XY = locXY == -1 ? v_double4_0 : cumSum[locXY];
            double4_t XZ = locXZ == -1 ? v_double4_0 : cumSum[locXZ];
            double4_t XYZW = cumSum[locXYZW];
            double4_t whole = cumSum[locWhole];

            double4_t sumInside = XYZW - XY - XZ + X;
            double4_t sumOutside = whole - sumInside;

            innerEach = 1.0 / v_nInside * sumInside;
            outerEach = 1.0 / v_nOutside * sumOutside;

            double4_t tmp = v_nInside * innerEach * innerEach - 2 * innerEach * sumInside +
                            v_nOutside * outerEach * outerEach - 2 * outerEach * sumOutside;

            double cost = sum_v(tmp);

            if (cost < bestCostEach)
            {
                bestCostEach = cost;
                bestResultEach.x0 = l;
                bestResultEach.y0 = t;
                bestResultEach.x1 = r;
                bestResultEach.y1 = b;

                bestResultEach.inner[0] = innerEach[0];
                bestResultEach.inner[1] = innerEach[1];
                bestResultEach.inner[2] = innerEach[2];
                bestResultEach.outer[0] = outerEach[0];
                bestResultEach.outer[1] = outerEach[1];
                bestResultEach.outer[2] = outerEach[2];

                // printf("size: %d -> (t, l, b, r) = (%d, %d, %d, %d): inner = [%.4f, %.4f, %.4f], outer = [%.4f, %.4f, %.4f] -> cost = %.4f\n", size, t, l, b + 1, r + 1, innerEach[0], innerEach[1], innerEach[2], outerEach[0], outerEach[1], outerEach[2], bestCostEach);
            }
        }

        costs[size] = bestCostEach;
        bestResult[size] = bestResultEach;
    }

    t2 = high_resolution_clock::now();
    ms_double = t2 - t1;
    printf("2. Running time: %8.4f\n", ms_double.count());

    for (int size = 1; size < NYX_pad; ++size)
    {
        int t = bestResult[size].y0;
        int l = bestResult[size].x0;
        int b = bestResult[size].y1 + 1;
        int r = bestResult[size].x1 + 1;
        float inner0 = bestResult[size].inner[0];
        float inner1 = bestResult[size].inner[1];
        float inner2 = bestResult[size].inner[2];
        float outer0 = bestResult[size].outer[0];
        float outer1 = bestResult[size].outer[1];
        float outer2 = bestResult[size].outer[2];
        // printf("(t, l, b, r) = (%d, %d, %d, %d): inner = [%.4f, %.4f, %.4f], outer = [%.4f, %.4f, %.4f] -> cost = %.4f\n", t, l, b + 1, r + 1, inner0, inner1, inner2, outer0, outer1, outer2, costs[size]);
    }

    // 3. Select best for each size

    t1 = high_resolution_clock::now();

    vector<int> bestIndices = vector<int>(N_VECS, 0);
    vector<double> bestCosts = vector<double>(N_VECS, inf);
    // TODO: HoangLe [May-25]: Padding this

#pragma omp parallel for schedule(dynamic, 1)
    for (int k = 0; k < N_VECS; ++k)
    {
        for (int size = n * k; size < n * (k + 1); ++size)
        {

            // NOTE: HoangLe [May-25]: Can apply instruction-level optimization here
            if (costs[size] < bestCosts[k])
            {
                bestCosts[k] = costs[size];
                bestIndices[k] = size;
            }
        }
        // printf("(t, l, b, r) = (%d, %d, %d, %d): inner = [%.4f, %.4f, %.4f], outer = [%.4f, %.4f, %.4f] -> cost = %.4f\n", t, l, b + 1, r + 1, inner[0], inner[1], inner[2], outer[0], outer[1], outer[2], cost);
    }

    t2 = high_resolution_clock::now();
    ms_double = t2 - t1;
    printf("3. Running time: %8.4f\n", ms_double.count());

    // 4. Choose best from best

    t1 = high_resolution_clock::now();

    double bestCost = inf;
    int bestIdx = 0;
    for (int k = 0; k < N_VECS; ++k)
    {
        if (bestCosts[k] < bestCost)
        {
            bestCost = bestCosts[k];
            bestIdx = bestIndices[k];
        }
    }

    result.y0 = bestResult[bestIdx].y0;
    result.x0 = bestResult[bestIdx].x0;
    result.y1 = bestResult[bestIdx].y1 + 1;
    result.x1 = bestResult[bestIdx].x1 + 1;

    result.inner[0] = bestResult[bestIdx].inner[0];
    result.inner[1] = bestResult[bestIdx].inner[1];
    result.inner[2] = bestResult[bestIdx].inner[2];
    result.outer[0] = bestResult[bestIdx].outer[0];
    result.outer[1] = bestResult[bestIdx].outer[1];
    result.outer[2] = bestResult[bestIdx].outer[2];

    t2 = high_resolution_clock::now();
    ms_double = t2 - t1;
    printf("4. Running time: %8.4f\n", ms_double.count());

    // Testing
    // Testing::printColor<vector<double>>(cumSum, nx, ny);

    return result;
}

int main(int argc, char const *argv[])
{
    // 1. Generate or read data from file
    int ny = 200, nx = 200;
    srand(0);
    Input input = {.ny = ny, .nx = nx, .data = DataGen::genRandArr(nx, ny)};

    // printColor(data, nx, ny);
    // Testing::printColor<float>(data, nx, ny);

    // string filename = FILE2;
    // Input input = DataGen::readFile(filename);
    // if (input.data == nullptr)
    // {
    //     cerr << "Cannot read data file: " << filename << endl;
    //     return 1;
    // }

    // DataGen::printInput(input);

    printf("\n===================\n\n");

    Result result = segment(input.ny, input.nx, input.data);

    printf("\n===================\n\n");

    Testing::printResult(result);

    delete[] input.data;

    return 0;
}

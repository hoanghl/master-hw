#include <cassert>
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

static inline int getNumCells(int left, int top, int right, int bot)
{
    int nRowsInside = bot - top + 1;
    int nColsInside = right - left + 1;

    return nRowsInside * nColsInside;
}

Result segment(int ny, int nx, const float *data)
{
    const int NYX = ny * nx;
    constexpr int N_VECS = 8;

    int n = (NYX + N_VECS - 1) / N_VECS;
    int NYX_pad = n * N_VECS;

    Result result{0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}};

    // 1. Pre-compute

    vector<double4_t> cumSum = vector<double4_t>(nx * ny, v_double4_0);
    vector<double> cumSum3 = vector<double>(nx * ny, 0);
    vector<double> cumSum3_2 = vector<double>(nx * ny, 0);
    vector<int> locationsX = vector<int>(ny * nx);
    vector<int> locationsXY = vector<int>(ny * nx);
    vector<int> locationsXZ = vector<int>(ny * nx);

    auto t1 = high_resolution_clock::now();

    cumSum[0][0] = data[0];
    cumSum[0][1] = data[1];
    cumSum[0][2] = data[2];
    cumSum3[0] = data[0] + data[1] + data[2];
    cumSum3_2[0] = -2 * (data[0] * data[1] + data[0] * data[2] + data[1] * data[2]);

    for (int br = 1; br < NYX; ++br)
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
        cumSum3[locCurrent] = sum_v(XY + XZ - X + current);
        cumSum3_2[locCurrent] = -2 * (cumSum[locCurrent][0] * cumSum[locCurrent][1] + cumSum[locCurrent][0] * cumSum[locCurrent][2] + cumSum[locCurrent][1] * cumSum[locCurrent][2]);
    }

#pragma omp parallel for collapse(2)
    for (int y = 0; y < ny; ++y)
        for (int x = 0; x < nx; ++x)
        {
            int idx = y * nx + x;

            locationsX[idx] = x == 0 || y == 0 ? -1 : (x - 1) + nx * (y - 1);
            locationsXY[idx] = y == 0 ? -1 : x + nx * (y - 1);
            locationsXZ[idx] = x == 0 ? -1 : (x - 1) + nx * y;
        }

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    printf("1. Running time: %8.4f\n", ms_double.count());

    // 2. Main loop

    t1 = high_resolution_clock::now();

    vector<double> costs = vector<double>(NYX_pad, inf);
    vector<TLBR> bestTLBR = vector<TLBR>(NYX_pad);

    int locWhole = nx * ny - 1;
    double4_t whole = cumSum[locWhole];

#pragma omp parallel for schedule(dynamic, ny) collapse(2)
    for (int height = 1; height <= ny; ++height)
        for (int width = 1; width <= nx; ++width)
        {

            int nInside = height * width;
            int nOutside = NYX - nInside;

            double bestCostEach = inf;
            int bestY0, bestX0, bestY1, bestX1;
            double cost;

            if (nOutside > 0)
            {
                for (int t = 0; t <= ny - height; ++t)
                    for (int l = 0; l <= nx - width; ++l)
                    {
                        // Determine b, r
                        int b = t + height - 1, r = l + width - 1;

                        // Calculate components for include-exclude principle
                        int locX = locationsX[t * nx + l];
                        int locXY = locationsXY[t * nx + r];
                        int locXZ = locationsXZ[b * nx + l];
                        int locXYZW = r + nx * b;

                        // double4_t X = locX == -1 ? v_double4_0 : cumSum[locX];
                        // double4_t XY = locXY == -1 ? v_double4_0 : cumSum[locXY];
                        // double4_t XZ = locXZ == -1 ? v_double4_0 : cumSum[locXZ];
                        // double4_t XYZW = cumSum[locXYZW];

                        // double4_t sumInside = XYZW - XY - XZ + X;
                        // double4_t sumOutside = whole - sumInside;

                        // double4_t v_nInside = _mm256_set1_pd(-1.0 / nInside);
                        // double4_t v_nOutside = _mm256_set1_pd(-1.0 / nOutside);

                        // cost = sum_v(
                        //     _mm256_add_pd(
                        //         _mm256_mul_pd(_mm256_mul_pd(v_nInside, sumInside), sumInside),
                        //         _mm256_mul_pd(_mm256_mul_pd(v_nOutside, sumOutside), sumOutside)));

                        double4_t X = locX == -1 ? v_double4_0 : cumSum[locX];
                        double4_t XY = locXY == -1 ? v_double4_0 : cumSum[locXY];
                        double4_t XZ = locXZ == -1 ? v_double4_0 : cumSum[locXZ];
                        double4_t XYZW = cumSum[locXYZW];

                        double4_t sumInside = XYZW - XY - XZ + X;
                        double4_t sumOutside = whole - sumInside;

                        double4_t v_nInside = _mm256_set1_pd(-1.0 / nInside);
                        double4_t v_nOutside = _mm256_set1_pd(-1.0 / nOutside);

                        cost = sum_v(
                            _mm256_add_pd(
                                _mm256_mul_pd(_mm256_mul_pd(v_nInside, sumInside), sumInside),
                                _mm256_mul_pd(_mm256_mul_pd(v_nOutside, sumOutside), sumOutside)));
                                
                        // cost = -1.0 / nInside * sum_v(_mm256_mul_pd(sumInside, sumInside))
                        // -1.0 / nOutside * sum_v(_mm256_mul_pd(sumOutside, sumOutside));

                        if (cost < bestCostEach)
                        {
                            bestCostEach = cost;

                            bestX0 = l;
                            bestY0 = t;
                            bestX1 = r;
                            bestY1 = b;
                        }
                    }
            }

            int idx = (height - 1) * nx + (width - 1);
            costs[idx] = bestCostEach;
            bestTLBR[idx] = std::make_tuple(bestY0, bestX0, bestY1, bestX1);

            // printf("idx = %3d: costs[%2d] = %.4f, (t, l, b, r) = (%2d, %2d, %2d, %2d)\n", idx, idx, costs[idx], bestY0, bestX0, bestY1, bestX1);
        }

    t2 = high_resolution_clock::now();
    ms_double = t2 - t1;
    printf("2. Running time: %8.4f\n", ms_double.count());

    // for (int i = 0; i < NYX_pad; ++i)
    // {
    //     int t = get<0>(bestTLBR[i]);
    //     int l = get<1>(bestTLBR[i]);
    //     int b = get<2>(bestTLBR[i]);
    //     int r = get<3>(bestTLBR[i]);

    //     printf("i = %3d: costs[%2d] = %.4f, (t, l, b, r) = (%2d, %2d, %2d, %2d)\n", i, i, costs[i], t, l, b, r);
    // }

    // 3. Select best for each size

    t1 = high_resolution_clock::now();

    vector<int> bestIndices = vector<int>(N_VECS, 0);
    vector<double> bestCosts = vector<double>(N_VECS, inf);

#pragma omp parallel for schedule(dynamic, 1)
    for (int k = 0; k < N_VECS; ++k)
    {
        for (int size = n * k; size < n * (k + 1); ++size)
        {
            if (costs[size] < bestCosts[k])
            {
                bestCosts[k] = costs[size];
                bestIndices[k] = size;
            }
        }
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

    // Calculate inner and outer

    int bestT = get<0>(bestTLBR[bestIdx]);
    int bestL = get<1>(bestTLBR[bestIdx]);
    int bestB = get<2>(bestTLBR[bestIdx]);
    int bestR = get<3>(bestTLBR[bestIdx]);

    int locX = locationsX[bestT * nx + bestL];
    int locXY = locationsXY[bestT * nx + bestR];
    int locXZ = locationsXZ[bestB * nx + bestL];
    int locXYZW = bestR + nx * bestB;

    double4_t X = locX == -1 ? v_double4_0 : cumSum[locX];
    double4_t XY = locXY == -1 ? v_double4_0 : cumSum[locXY];
    double4_t XZ = locXZ == -1 ? v_double4_0 : cumSum[locXZ];
    double4_t XYZW = cumSum[locXYZW];
    double4_t sumInside = XYZW - XY - XZ + X;
    double4_t sumOutside = whole - sumInside;

    int nInside = (bestB - bestT + 1) * (bestR - bestL + 1);
    int nOutside = NYX - nInside;

    double4_t v_inner = 1.0 / nInside * sumInside;
    double4_t v_outer = 1.0 / nOutside * sumOutside;

    result.y0 = bestT;
    result.x0 = bestL;
    result.y1 = bestB + 1;
    result.x1 = bestR + 1;

    result.inner[0] = static_cast<float>(v_inner[0]);
    result.inner[1] = static_cast<float>(v_inner[1]);
    result.inner[2] = static_cast<float>(v_inner[2]);
    result.outer[0] = static_cast<float>(v_outer[0]);
    result.outer[1] = static_cast<float>(v_outer[1]);
    result.outer[2] = static_cast<float>(v_outer[2]);

    t2 = high_resolution_clock::now();
    ms_double = t2 - t1;
    printf("4. Running time: %8.4f\n", ms_double.count());

    return result;
}

int main(int argc, char const *argv[])
{
    // 1. Generate or read data from file
    // srand(0);
    int ny = 400, nx = 400;
    Input input = {.ny = ny, .nx = nx, .data = DataGen::genRandArr(nx, ny)};

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

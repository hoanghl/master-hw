#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <tuple>
#include <vector>

#include "exp.hpp"

using namespace std;
using namespace std::chrono;

Result segment(int ny, int nx, const float *data)
{
    // NOTE: HoangLe [May-24]: In my code, bottom-r is in included mode whereas the requirement of the output is excluded mode

    const int NYX = ny * nx;
    constexpr int N_VECS = 8;

    int n = (NYX * NYX + N_VECS - 1) / N_VECS;
    const int NYX2 = n * N_VECS;

    Result result{0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}};

    // 1. Pre-compute

    // 1.1. Define pre-computed vectors
    vector<double> cumSum = vector<double>(N_COLOR_CHAN * nx * ny, 0.);
    vector<double> cumSumSq = vector<double>(N_COLOR_CHAN * nx * ny, 0.);
    cumSum[0] = data[0];
    cumSum[1] = data[1];
    cumSum[2] = data[2];

    // 1.2. Pre-compute
    auto t1 = high_resolution_clock::now();

    // NOTE: HoangLe [May-21]: Can use multithread, z-order here

    for (int br = 1; br < NYX; ++br) // Start at 2 because position [0, 0] is pre-assigned above
    {
        int r, b;
        Utils::xy2each(br, nx, r, b);
        int t = b - 1, l = r - 1;

        for (int c = 0; c < N_COLOR_CHAN; ++c)
        {
            int locX = t < 0 || l < 0 ? -1 : c + 3 * (l + nx * t);
            int locXY = t < 0 ? -1 : c + 3 * (r + nx * t);
            int locXZ = l < 0 ? -1 : c + 3 * (l + nx * b);
            int locCurrent = c + 3 * (r + nx * b);

            double X = locX == -1 ? 0 : cumSum[locX];
            double XY = locXY == -1 ? 0 : cumSum[locXY];
            double XZ = locXZ == -1 ? 0 : cumSum[locXZ];
            double current = data[locCurrent];

            cumSum[locCurrent] = XY + XZ - X + current;
        }
    }

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    printf("1. Running time: %8.4f\n", ms_double.count());

    // 2. Assign inner, outer, cost

    t1 = high_resolution_clock::now();

    vector<COLOR> inner = vector<COLOR>(NYX2, make_tuple(inf, inf, inf));
    vector<COLOR> outer = vector<COLOR>(NYX2, make_tuple(inf, inf, inf));
    vector<double> costs = vector<double>(NYX2, inf);

    // TODO: HoangLe [May-25]: Fill 3 above vectors with heavy optimizations

#pragma omp parallel for schedule(dynamic, ny)
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

        if (nOutside == 0)
            continue;

        double innerEach[N_COLOR_CHAN] = {};
        double outerEach[N_COLOR_CHAN] = {};
        costs[tlbr] = 0;

        for (int c = 0; c < N_COLOR_CHAN; ++c)
        {
            // Calculate components for include-exclude principle

            int locX = l == 0 || t == 0 ? -1 : c + 3 * ((l - 1) + nx * (t - 1));
            int locXY = t == 0 ? -1 : c + 3 * (r + nx * (t - 1));
            int locXZ = l == 0 ? -1 : c + 3 * ((l - 1) + nx * b);
            int locXYZW = c + 3 * (r + nx * b);
            int locWhole = c + 3 * (nx * ny - 1);

            double X = locX == -1 ? 0 : cumSum[locX];
            double XY = locXY == -1 ? 0 : cumSum[locXY];
            double XZ = locXZ == -1 ? 0 : cumSum[locXZ];
            double XYZW = cumSum[locXYZW];
            double whole = cumSum[locWhole];

            double Xsq = locX == -1 ? 0 : cumSumSq[locX];
            double XYsq = locXY == -1 ? 0 : cumSumSq[locXY];
            double XZsq = locXZ == -1 ? 0 : cumSumSq[locXZ];
            double XYZWsq = cumSumSq[locXYZW];
            double wholesq = cumSumSq[locWhole];

            double sumInside = XYZW - XY - XZ + X;
            double sumOutside = whole - sumInside;

            // Calculate inner and outer
            innerEach[c] = 1.0 / nInside * sumInside;
            outerEach[c] = 1.0 / nOutside * sumOutside;

            // Calculate cost
            costs[tlbr] += nInside * pow(innerEach[c], 2) - 2 * innerEach[c] * sumInside;
            costs[tlbr] += nOutside * pow(outerEach[c], 2) - 2 * outerEach[c] * sumOutside;
        }

        inner[tlbr] = make_tuple(innerEach[0], innerEach[1], innerEach[2]);
        outer[tlbr] = make_tuple(outerEach[0], outerEach[1], outerEach[2]);
    }

    t2 = high_resolution_clock::now();
    ms_double = t2 - t1;
    printf("2. Running time: %8.4f\n", ms_double.count());

    // for (int tlbr = 0; tlbr < NYX2; ++tlbr)
    // {
    //     int tl, br, t, l, b, r;
    //     tl = tlbr / NYX;
    //     br = tlbr % NYX;
    //     t = tl / nx;
    //     l = tl % nx;
    //     b = br / nx;
    //     r = br % nx;

    //     if (b < t || r < l || t >= ny || b >= ny || l >= nx || r >= nx)
    //         continue;

    //     double inner0 = get<0>(inner[tlbr]);
    //     double inner1 = get<1>(inner[tlbr]);
    //     double inner2 = get<2>(inner[tlbr]);
    //     double outer0 = get<0>(outer[tlbr]);
    //     double outer1 = get<1>(outer[tlbr]);
    //     double outer2 = get<2>(outer[tlbr]);
    //     printf("(t, l, b, r) = (%d, %d, %d, %d): inner = [%.4f, %.4f, %.4f], outer = [%.4f, %.4f, %.4f] -> cost = %.4f\n", t, l, b + 1, r + 1, inner0, inner1, inner2, outer0, outer1, outer2, costs[tlbr]);
    // }

    // 3. Select best for each size

    t1 = high_resolution_clock::now();

    vector<int> bestIndices = vector<int>(N_VECS, 0);
    vector<double> bestCosts = vector<double>(N_VECS, inf);
    // TODO: HoangLe [May-25]: Padding this

#pragma omp parallel for schedule(dynamic, 1)
    for (int k = 0; k < N_VECS; ++k)
    {
        for (int tlbr = n * k; tlbr < n * (k + 1); ++tlbr)
        {

            // NOTE: HoangLe [May-25]: Can apply instruction-level optimization here
            if (costs[tlbr] < bestCosts[k])
            {
                bestCosts[k] = costs[tlbr];
                bestIndices[k] = tlbr;
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

    int tl, br;
    tl = bestIdx / NYX;
    br = bestIdx % NYX;
    result.y0 = tl / nx;
    result.x0 = tl % nx;
    result.y1 = br / nx + 1;
    result.x1 = br % nx + 1;

    result.inner[0] = static_cast<float>(get<0>(inner[bestIdx]));
    result.inner[1] = static_cast<float>(get<1>(inner[bestIdx]));
    result.inner[2] = static_cast<float>(get<2>(inner[bestIdx]));
    result.outer[0] = static_cast<float>(get<0>(outer[bestIdx]));
    result.outer[1] = static_cast<float>(get<1>(outer[bestIdx]));
    result.outer[2] = static_cast<float>(get<2>(outer[bestIdx]));

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

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
    // NOTE: HoangLe [May-24]: In my code, bottom-right is in included mode whereas the requirement of the output is excluded mode

    Result result{0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}};

    // 1. Pre-compute

    // 1.1. Define pre-computed vectors
    vector<double> cumSum = vector<double>(N_COLOR_CHAN * (nx + 1) * (ny + 1), 0.);
    vector<double> cumSumSq = vector<double>(N_COLOR_CHAN * (nx + 1) * (ny + 1), 0.);
    cumSum[0] = data[0];
    cumSum[1] = data[1];
    cumSum[2] = data[2];
    cumSumSq[0] = pow(data[0], 2);
    cumSumSq[1] = pow(data[1], 2);
    cumSumSq[2] = pow(data[2], 2);

    // 1.2. Pre-compute
    auto t1 = high_resolution_clock::now();

    // NOTE: HoangLe [May-21]: Can use multithread, z-order here

    // #pragma omp parallel for schedule(static, 1)
    for (int br = 1; br < ny * nx; ++br) // Start at 2 because position [0, 0] is pre-assigned above
    {
        int right, bot;
        Utils::xy2each(br, nx, right, bot);
        int top = bot - 1, left = right - 1;

        for (int c = 0; c < N_COLOR_CHAN; ++c)
        {
            int locX = top < 0 || left < 0 ? -1 : c + 3 * (left + nx * top);
            int locXY = top < 0 ? -1 : c + 3 * (right + nx * top);
            int locXZ = left < 0 ? -1 : c + 3 * (left + nx * bot);
            int locCurrent = c + 3 * (right + nx * bot);

            double X = locX == -1 ? 0 : cumSum[locX];
            double XY = locXY == -1 ? 0 : cumSum[locXY];
            double XZ = locXZ == -1 ? 0 : cumSum[locXZ];
            double Xsq = locX == -1 ? 0 : cumSumSq[locX];
            double XYsq = locXY == -1 ? 0 : cumSumSq[locXY];
            double XZsq = locXZ == -1 ? 0 : cumSumSq[locXZ];

            double current = data[locCurrent];
            double currentSq = pow(data[locCurrent], 2);

            cumSum[locCurrent] = XY + XZ - X + current;
            cumSumSq[locCurrent] = XYsq + XZsq - Xsq + currentSq;
        }
    }

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    printf("1. Running time: %8.4f\n", ms_double.count());

    // 2. Start

    t1 = high_resolution_clock::now();

    vector<Result> bestEachSize = vector<Result>(ny * nx);
    vector<double> costEachSize = vector<double>(ny * nx, inf);

// NOTE: HoangLe [May-21]: Can use multithread here
#pragma omp parallel for schedule(static, 1)
    for (int size = 1; size < ny * nx; ++size)
    {
        Result best{0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}};
        double bestCost = inf;

        // printf("size = %d\n", size);

        for (int tl = 0; tl < ny * nx; ++tl)
        {
            int top = 0, left = 0;
            Utils::xy2each(tl, nx, left, top);

            for (int bot = top; bot < min(top + size + 1, ny); ++bot)
            {
                for (int right = left; right < min(left + size + 1, nx); ++right)
                {

                    if (Utils::getNumCells(left, top, right, bot) != size)
                    {
                        continue;
                    }

                    vector<double> inner = vector<double>(3, 0.);
                    vector<double> outer = vector<double>(3, 0.);
                    double cost = 0;

                    int nInside = Utils::getNumCells(left, top, right, bot);
                    int nOutside = Utils::getNumCellsOutside(left, top, right, bot, nx, ny);

                    for (int c = 0; c < N_COLOR_CHAN; ++c)
                    {
                        // Calculate components for include-exclude principle

                        int locX = left == 0 || top == 0 ? -1 : c + 3 * ((left - 1) + nx * (top - 1));
                        int locXY = top == 0 ? -1 : c + 3 * (right + nx * (top - 1));
                        int locXZ = left == 0 ? -1 : c + 3 * ((left - 1) + nx * bot);
                        int locXYZW = c + 3 * (right + nx * bot);
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
                        double sumInsideSq = XYZWsq - XYsq - XZsq + Xsq;
                        double sumOutsideSq = wholesq - sumInsideSq;

                        // Calculate inner and outer
                        inner[c] = 1.0 / nInside * sumInside;
                        outer[c] = 1.0 / nOutside * sumOutside;

                        // Calculate cost
                        cost += nInside * pow(inner[c], 2) - 2 * inner[c] * sumInside + sumInsideSq;
                        cost += nOutside * pow(outer[c], 2) - 2 * outer[c] * sumOutside + sumOutsideSq;
                    }

                    // printf("(t, l, b, r) = (%d, %d, %d, %d): inner = [%.4f, %.4f, %.4f], outer = [%.4f, %.4f, %.4f] -> cost = %.4f\n", top, left, bot + 1, right + 1, inner[0], inner[1], inner[2], outer[0], outer[1], outer[2], cost);

                    if (cost < bestCost)
                    {
                        bestCost = cost;

                        best.x0 = left;
                        best.y0 = top;
                        best.x1 = right + 1;
                        best.y1 = bot + 1;
                        best.inner[0] = inner[0];
                        best.inner[1] = inner[1];
                        best.inner[2] = inner[2];
                        best.outer[0] = outer[0];
                        best.outer[1] = outer[1];
                        best.outer[2] = outer[2];
                    }
                }
            }
        }

        bestEachSize[size] = best;
        costEachSize[size] = bestCost;
    }

    t2 = high_resolution_clock::now();
    ms_double = t2 - t1;
    printf("2. Running time: %8.4f\n", ms_double.count());

    // 3. Choose best tl, br

    t1 = high_resolution_clock::now();

    // NOTE: HoangLe [May-21]: Can use multithread here
    double bestCost = inf;
    double bestIdx = 0;
    for (int i = 1; i < ny * nx; ++i)
    {
        if (costEachSize[i] < bestCost)
        {
            bestCost = costEachSize[i];
            bestIdx = i;
        }
    }

    result.x0 = bestEachSize[bestIdx].x0;
    result.y0 = bestEachSize[bestIdx].y0;
    result.x1 = bestEachSize[bestIdx].x1;
    result.y1 = bestEachSize[bestIdx].y1;
    result.inner[0] = static_cast<float>(bestEachSize[bestIdx].inner[0]);
    result.inner[1] = static_cast<float>(bestEachSize[bestIdx].inner[1]);
    result.inner[2] = static_cast<float>(bestEachSize[bestIdx].inner[2]);
    result.outer[0] = static_cast<float>(bestEachSize[bestIdx].outer[0]);
    result.outer[1] = static_cast<float>(bestEachSize[bestIdx].outer[1]);
    result.outer[2] = static_cast<float>(bestEachSize[bestIdx].outer[2]);

    t2 = high_resolution_clock::now();
    ms_double = t2 - t1;
    printf("3. Running time: %8.4f\n", ms_double.count());

    // Testing
    // Testing::printColor<vector<double>>(cumSum, nx, ny);

    return result;
}

int main(int argc, char const *argv[])
{
    // 1. Generate or read data from file
    int ny = 100, nx = 100;
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

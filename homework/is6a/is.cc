// #include <cstdlib>
#include <limits>
#include <math.h>
#include <tuple>
#include <vector>

using namespace std;

typedef tuple<int, int, int, int> TLBR;

const float inf = numeric_limits<float>::infinity();

struct Result
{
    int y0;
    int x0;
    int y1;
    int x1;
    float outer[3];
    float inner[3];
};

/*
This is the function you need to implement. Quick reference:
- x coordinates: 0 <= x < nx
- y coordinates: 0 <= y < ny
- color components: 0 <= c < 3
- input: data[c + 3 * x + 3 * nx * y]
*/
Result segment(int ny, int nx, const float *data)
{
    const int NYX = ny * nx;
    constexpr int N_VECS = 8;

    int n = (NYX + N_VECS - 1) / N_VECS;
    int NYX_pad = n * N_VECS;

    Result result{0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}};

    // 1. Pre-compute

    vector<int> cumSum = vector<int>(nx * ny, 0);
    vector<int> locationsX = vector<int>(ny * nx);
    vector<int> locationsXY = vector<int>(ny * nx);
    vector<int> locationsXZ = vector<int>(ny * nx);

    cumSum[0] = static_cast<int>(data[0]);

    for (int br = 1; br < NYX; ++br)
    {
        int r = br % nx, b = br / nx;
        int t = b - 1, l = r - 1;

        // Determine regions' location
        int locX = t < 0 || l < 0 ? -1 : l + nx * t;
        int locXY = t < 0 ? -1 : r + nx * t;
        int locXZ = l < 0 ? -1 : l + nx * b;
        int locCurrent = r + nx * b;

        int X = locX == -1 ? 0 : cumSum[locX];
        int XY = locXY == -1 ? 0 : cumSum[locXY];
        int XZ = locXZ == -1 ? 0 : cumSum[locXZ];
        int current = static_cast<int>(data[3 * locCurrent]);

        cumSum[locCurrent] = XY + XZ - X + current;
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

    // 2. Assign inner, outer, cost

    vector<float> costs = vector<float>(NYX_pad, inf);
    vector<TLBR> bestTLBR = vector<TLBR>(NYX_pad);

    int locWhole = nx * ny - 1;
    int whole = cumSum[locWhole];

#pragma omp parallel for collapse(2) schedule(dynamic, ny)
    for (int height = 1; height <= ny; ++height)
        for (int width = 1; width <= nx; ++width)
        {
            int nInside = height * width;
            int nOutside = NYX - nInside;

            float bestCostEach = inf;
            int bestY0 = 0, bestX0 = 0, bestY1 = 0, bestX1 = 0;

            if (nOutside > 0)
            {
                for (int t = 0; t <= ny - height; ++t)
                {
                    for (int l = 0; l <= nx - width; ++l)
                    {
                        int b = t + height - 1;
                        int r = l + width - 1;

                        int locX = locationsX[t * nx + l];
                        int locXY = locationsXY[t * nx + r];
                        int locXZ = locationsXZ[b * nx + l];
                        int locXYZW = r + nx * b;

                        int X = locX == -1 ? 0 : cumSum[locX];
                        int XY = locXY == -1 ? 0 : cumSum[locXY];
                        int XZ = locXZ == -1 ? 0 : cumSum[locXZ];
                        int XYZW = cumSum[locXYZW];

                        int sumInside = XYZW - XY - XZ + X;
                        int sumOutside = whole - sumInside;

                        float cost = -1.0 / nInside * powf(sumInside, 2) - 1.0 / nOutside * powf(sumOutside, 2);

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
            }

            int idx = (height - 1) * nx + (width - 1);
            costs[idx] = bestCostEach;
            bestTLBR[idx] = std::make_tuple(bestY0, bestX0, bestY1, bestX1);
        }

    // 3. Select best for each size

    vector<int> bestIndices = vector<int>(N_VECS, 0);
    vector<float> bestCosts = vector<float>(N_VECS, inf);

#pragma omp parallel for schedule(dynamic, 1)
    for (int k = 0; k < N_VECS; ++k)
    {
        for (int tlbr = n * k; tlbr < n * (k + 1); ++tlbr)
        {
            if (costs[tlbr] < bestCosts[k])
            {
                bestCosts[k] = costs[tlbr];
                bestIndices[k] = tlbr;
            }
        }
    }

    // 4. Choose best from best
    float bestCost = inf;
    int bestIdx = 0;
    for (int k = 0; k < N_VECS; ++k)
    {
        if (bestCosts[k] < bestCost)
        {
            bestCost = bestCosts[k];
            bestIdx = bestIndices[k];
        }
    }

    int bestT = get<0>(bestTLBR[bestIdx]);
    int bestL = get<1>(bestTLBR[bestIdx]);
    int bestB = get<2>(bestTLBR[bestIdx]);
    int bestR = get<3>(bestTLBR[bestIdx]);

    int locX = bestL == 0 || bestT == 0 ? -1 : (bestL - 1) + nx * (bestT - 1);
    int locXY = bestT == 0 ? -1 : bestR + nx * (bestT - 1);
    int locXZ = bestL == 0 ? -1 : (bestL - 1) + nx * bestB;
    int locXYZW = bestR + nx * bestB;

    int X = locX == -1 ? 0 : cumSum[locX];
    int XY = locXY == -1 ? 0 : cumSum[locXY];
    int XZ = locXZ == -1 ? 0 : cumSum[locXZ];
    int XYZW = cumSum[locXYZW];
    int sumInside = XYZW - XY - XZ + X;
    int sumOutside = whole - sumInside;

    int nInside = (bestB - bestT + 1) * (bestR - bestL + 1);
    int nOutside = NYX - nInside;

    float inner = 1.0 / nInside * sumInside;
    float outer = 1.0 / nOutside * sumOutside;

    result.y0 = bestT;
    result.x0 = bestL;
    result.y1 = bestB + 1;
    result.x1 = bestR + 1;

    result.inner[0] = inner;
    result.inner[1] = inner;
    result.inner[2] = inner;
    result.outer[0] = outer;
    result.outer[1] = outer;
    result.outer[2] = outer;

    return result;
}

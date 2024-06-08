// #include <cstdlib>
#include <limits>
#include <math.h>
#include <tuple>
#include <vector>
#include <x86intrin.h>

using namespace std;

typedef double double4_t __attribute__((vector_size(4 * sizeof(double))));
typedef tuple<int, int, int, int> TLBR;

const double4_t v_double4_0 = {
    0.,
    0.,
    0.,
    0.,
};

static inline double sum_v(double4_t x) { return x[0] + x[1] + x[2]; }

const double inf = numeric_limits<double>::infinity();

static inline int getNumCells(int left, int top, int right, int bot)
{
    int nRowsInside = bot - top + 1;
    int nColsInside = right - left + 1;

    return nRowsInside * nColsInside;
}

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
    // NOTE: HoangLe [May-24]: In my code, bottom-right is in included mode whereas the requirement of the output is excluded mode

    const int NYX = ny * nx;
    constexpr int N_VECS = 8;

    int n = (NYX + N_VECS - 1) / N_VECS;
    int NYX_pad = n * N_VECS;

    Result result{0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}};

    // 1. Pre-compute

    vector<double4_t> cumSum = vector<double4_t>(nx * ny, v_double4_0);
    // vector<double> cumSum3 = vector<double>(nx * ny, 0);
    vector<int> locationsX = vector<int>(ny * nx);
    vector<int> locationsXY = vector<int>(ny * nx);
    vector<int> locationsXZ = vector<int>(ny * nx);

    cumSum[0][0] = data[0];
    cumSum[0][1] = data[1];
    cumSum[0][2] = data[2];
    // cumSum3[0] = data[0] + data[1] + data[2];

    for (int br = 1; br < NYX; ++br) // Start at 2 because position [0, 0] is pre-assigned above
    {
        int r = br % nx, b = br / nx;
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
        // cumSum3[locCurrent] = sum_v(XY + XZ - X + current);
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

    vector<double> costs = vector<double>(NYX_pad, inf);
    vector<TLBR> bestTLBR = vector<TLBR>(NYX_pad);

#pragma omp parallel for collapse(2) schedule(dynamic, ny)
    for (int height = 1; height <= ny; ++height)
        for (int width = 1; width <= nx; ++width)
        {
            int locWhole = nx * ny - 1;
            int nInside = height * width;
            int nOutside = NYX - nInside;
            // double whole3 = cumSum3[locWhole];

            double bestCostEach = inf;
            int bestY0 = 0, bestX0 = 0, bestY1 = 0, bestX1 = 0;
            double cost = 0.0;

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

                        double4_t X = locX == -1 ? v_double4_0 : cumSum[locX];
                        double4_t XY = locXY == -1 ? v_double4_0 : cumSum[locXY];
                        double4_t XZ = locXZ == -1 ? v_double4_0 : cumSum[locXZ];
                        double4_t XYZW = cumSum[locXYZW];
                        double4_t whole = cumSum[locWhole];

                        double4_t sumInside = XYZW - XY - XZ + X;
                        double4_t sumOutside = whole - sumInside;

                        double4_t v_nInside = _mm256_set1_pd(-1.0 / nInside);
                        double4_t v_nOutside = _mm256_set1_pd(-1.0 / nOutside);

                        cost = sum_v(
                            _mm256_add_pd(
                                _mm256_mul_pd(_mm256_mul_pd(v_nInside, sumInside), sumInside),
                                _mm256_mul_pd(_mm256_mul_pd(v_nOutside, sumOutside), sumOutside)));

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
    vector<double> bestCosts = vector<double>(N_VECS, inf);

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

    int bestT = get<0>(bestTLBR[bestIdx]);
    int bestL = get<1>(bestTLBR[bestIdx]);
    int bestB = get<2>(bestTLBR[bestIdx]);
    int bestR = get<3>(bestTLBR[bestIdx]);

    int locX = bestL == 0 || bestT == 0 ? -1 : (bestL - 1) + nx * (bestT - 1);
    int locXY = bestT == 0 ? -1 : bestR + nx * (bestT - 1);
    int locXZ = bestL == 0 ? -1 : (bestL - 1) + nx * bestB;
    int locXYZW = bestR + nx * bestB;

    int locWhole = nx * ny - 1;
    double4_t X = locX == -1 ? v_double4_0 : cumSum[locX];
    double4_t XY = locXY == -1 ? v_double4_0 : cumSum[locXY];
    double4_t XZ = locXZ == -1 ? v_double4_0 : cumSum[locXZ];
    double4_t XYZW = cumSum[locXYZW];
    double4_t sumInside = XYZW - XY - XZ + X;
    double4_t whole = cumSum[locWhole];
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

    return result;
}

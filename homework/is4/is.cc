#include <cstdlib>
#include <tuple>
#include <vector>

#include <x86intrin.h>

#include <math.h>

using namespace std;

typedef tuple<double, double, double> COLOR;

typedef double double4_t __attribute__((vector_size(4 * sizeof(double))));

const double4_t v_double4_0 = {
    0.,
    0.,
    0.,
    0.,
};

static inline double sum_v(double4_t x) { return x[0] + x[1] + x[2]; }

const double inf = numeric_limits<double>::infinity();

static inline int getNumCells(int x, int y, bool excludeLast = true)
{
    int nCols = excludeLast ? x : x + 1;
    int nRows = excludeLast ? y : y + 1;

    return nRows * nCols;
}

static inline int getNumCells(int left, int top, int right, int bot)
{
    int nRowsInside = bot - top + 1;
    int nColsInside = right - left + 1;

    return nRowsInside * nColsInside;
}

static inline int getNumCellsOutside(int left, int top, int right, int bot, int nx, int ny)
{
    return getNumCells(0, 0, nx - 1, ny - 1) - getNumCells(left, top, right, bot);
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
    const int NYX_pad = n * N_VECS;
    const int NYX2 = n * N_VECS;

    Result result{0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}};

    // 1. Pre-compute

    // 1.1. Define pre-computed vectors
    vector<double4_t> cumSum = vector<double4_t>(nx * ny, v_double4_0);
    // cumSum[0] = {data[0], data[1], data[2], 0.};
    cumSum[0][0] = data[0];
    cumSum[0][1] = data[1];
    cumSum[0][2] = data[2];

    // 1.2. Pre-compute

    // NOTE: HoangLe [May-21]: Can use multithread, z-order here

    // #pragma omp parallel for schedule(static, 1)
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
    }
    // 2. Assign inner, outer, cost

    vector<COLOR> inner = vector<COLOR>(NYX_pad, make_tuple(inf, inf, inf));
    vector<COLOR> outer = vector<COLOR>(NYX_pad, make_tuple(inf, inf, inf));
    vector<double> costs = vector<double>(NYX_pad, inf);

#pragma omp parallel for schedule(dynamic, ny)
    for (int size = 1; size < NYX_pad; ++size)
    {
        for (int tlbr = 0; tlbr < NYX2; ++tlbr)
        {
            int tl, br, t, l, b, r;
            tl = tlbr / NYX_pad;
            br = tlbr % NYX_pad;
            t = tl / nx;
            l = tl % nx;
            b = br / nx;
            r = br % nx;

            if (b < t || r < l || t >= ny || b >= ny || l >= nx || r >= nx)
                continue;

            int nInside = getNumCells(l, t, r, b);
            int nOutside = getNumCellsOutside(l, t, r, b, nx, ny);

            if (nOutside == 0 || nInside != size)
                continue;

            double4_t v_nInside = _mm256_set1_pd(nInside);
            double4_t v_nOutside = _mm256_set1_pd(nOutside);

            double4_t innerEach = v_double4_0;
            double4_t outerEach = v_double4_0;
            costs[tlbr] = 0;

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

            if (cost < costs[size])
            {
                costs[size] = cost;
                inner[size] = make_tuple(innerEach[0], innerEach[1], innerEach[2]);
                outer[size] = make_tuple(outerEach[0], outerEach[1], outerEach[2]);
            }
        }
    }

    // 3. Select best for each size

    vector<int> bestIndices = vector<int>(N_VECS, 0);
    vector<double> bestCosts = vector<double>(N_VECS, inf);

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

    // 4. Choose best from best

    // NOTE: HoangLe [May-21]: Can use multithread here
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

    return result;
}

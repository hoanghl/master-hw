// #include <cstdlib>
#include <vector>

#include <x86intrin.h>

#include <math.h>

using namespace std;

typedef double double4_t __attribute__((vector_size(4 * sizeof(double))));

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

    int n = (NYX * NYX + N_VECS - 1) / N_VECS;
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

    vector<double4_t> inner = vector<double4_t>(NYX2);
    vector<double4_t> outer = vector<double4_t>(NYX2);
    vector<double> costs = vector<double>(NYX2, inf);

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

        int nInside = getNumCells(l, t, r, b);
        int nOutside = getNumCellsOutside(l, t, r, b, nx, ny);

        if (nOutside == 0)
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

        costs[tlbr] = sum_v(tmp);
        inner[tlbr] = innerEach;
        outer[tlbr] = outerEach;
    }

    // 3. Select best for each size

    vector<int> bestIndices = vector<int>(N_VECS, 0);
    vector<double> bestCosts = vector<double>(N_VECS, inf);

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

    result.inner[0] = static_cast<float>(inner[bestIdx][0]);
    result.inner[1] = static_cast<float>(inner[bestIdx][1]);
    result.inner[2] = static_cast<float>(inner[bestIdx][2]);
    result.outer[0] = static_cast<float>(outer[bestIdx][0]);
    result.outer[1] = static_cast<float>(outer[bestIdx][1]);
    result.outer[2] = static_cast<float>(outer[bestIdx][2]);

    return result;
}

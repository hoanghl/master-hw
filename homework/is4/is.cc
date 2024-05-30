#include <cstdlib>
#include <tuple>
#include <vector>

#include <math.h>

using namespace std;

typedef tuple<int, int, int, int> TLBR;
typedef tuple<double, double, double> COLOR;

constexpr int N_COLOR_CHAN = 3;
const double inf = numeric_limits<double>::infinity();

void xy2each(int xy, int nx, int &x, int &y)
{
    x = xy % nx;
    y = xy / nx;
}

int getNumCells(int x, int y, bool excludeLast = true)
{
    int nCols = excludeLast ? x : x + 1;
    int nRows = excludeLast ? y : y + 1;

    return nRows * nCols;
}

int getNumCells(int left, int top, int right, int bot)
{
    int nRowsInside = bot - top + 1;
    int nColsInside = right - left + 1;

    return nRowsInside * nColsInside;
    ;
}

int getNumCellsOutside(int left, int top, int right, int bot, int nx, int ny)
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
    vector<double> cumSum = vector<double>(N_COLOR_CHAN * (nx + 1) * (ny + 1), 0.);
    // vector<double> cumSumSq = vector<double>(N_COLOR_CHAN * (nx + 1) * (ny + 1), 0.);
    cumSum[0] = data[0];
    cumSum[1] = data[1];
    cumSum[2] = data[2];
    // cumSumSq[0] = pow(data[0], 2);
    // cumSumSq[1] = pow(data[1], 2);
    // cumSumSq[2] = pow(data[2], 2);

    // 1.2. Pre-compute

    // NOTE: HoangLe [May-21]: Can use multithread, z-order here

    // #pragma omp parallel for schedule(static, 1)
    for (int br = 1; br < ny * nx; ++br) // Start at 2 because position [0, 0] is pre-assigned above
    {
        int right, bot;
        xy2each(br, nx, right, bot);
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
            // double Xsq = locX == -1 ? 0 : cumSumSq[locX];
            // double XYsq = locXY == -1 ? 0 : cumSumSq[locXY];
            // double XZsq = locXZ == -1 ? 0 : cumSumSq[locXZ];

            double current = data[locCurrent];
            // double currentSq = pow(data[locCurrent], 2);

            cumSum[locCurrent] = XY + XZ - X + current;
            // cumSumSq[locCurrent] = XYsq + XZsq - Xsq + currentSq;
        }
    }

    // 2. Assign inner, outer, cost

    vector<COLOR> inner = vector<COLOR>(NYX2, make_tuple(inf, inf, inf));
    vector<COLOR> outer = vector<COLOR>(NYX2, make_tuple(inf, inf, inf));
    vector<double> costs = vector<double>(NYX2, inf);

    // TODO: HoangLe [May-25]: Fill 3 above vectors with heavy optimizations

#pragma omp parallel for schedule(dynamic, 1000)
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

        // vector<double> innerEach = vector<double>(N_COLOR_CHAN, 0.);
        // vector<double> outerEach = vector<double>(N_COLOR_CHAN, 0.);
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

            // double Xsq = locX == -1 ? 0 : cumSumSq[locX];
            // double XYsq = locXY == -1 ? 0 : cumSumSq[locXY];
            // double XZsq = locXZ == -1 ? 0 : cumSumSq[locXZ];
            // double XYZWsq = cumSumSq[locXYZW];
            // double wholesq = cumSumSq[locWhole];

            double sumInside = XYZW - XY - XZ + X;
            double sumOutside = whole - sumInside;
            // double sumInsideSq = XYZWsq - XYsq - XZsq + Xsq;
            // double sumOutsideSq = wholesq - sumInsideSq;

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

    result.inner[0] = static_cast<float>(get<0>(inner[bestIdx]));
    result.inner[1] = static_cast<float>(get<1>(inner[bestIdx]));
    result.inner[2] = static_cast<float>(get<2>(inner[bestIdx]));
    result.outer[0] = static_cast<float>(get<0>(outer[bestIdx]));
    result.outer[1] = static_cast<float>(get<1>(outer[bestIdx]));
    result.outer[2] = static_cast<float>(get<2>(outer[bestIdx]));

    return result;
}

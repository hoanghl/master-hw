#include <cstdlib>
#include <tuple>
#include <vector>

using namespace std;

typedef tuple<int, int, int, int> TLBR;
constexpr int N_COLOR_CHAN = 3;
const double inf = numeric_limits<double>::infinity();

float genColor()
{
    return (rand() % 256) * 1.0;
}
float *genRandArr(int nx, int ny)
{
    float *out = new float[nx * ny * N_COLOR_CHAN];

    for (int y = 0; y < ny; ++y)
        for (int x = 0; x < nx; ++x)
            for (int c = 0; c < N_COLOR_CHAN; ++c)
            {
                out[c + 3 * x + 3 * nx * y] = genColor();
            }

    return out;
}

// == namespace: Utils =================================================

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

int getNumCells(int left, int top, int right, int bottom)
{
    int nRowsInside = bottom - top;
    int nColsInside = right - left;

    return nRowsInside * nColsInside;
    ;
}

int getNumCellsOutside(int left, int top, int right, int bottom, int nx, int ny)
{
    return getNumCells(nx, ny) - getNumCells(left, top, right, bottom);
}

bool isInside(int x, int y, int tl, int br, int nx)
{
    int left, top, right, bottom;
    xy2each(tl, nx, left, top);
    xy2each(br, nx, right, bottom);
    ++right;
    ++bottom;

    return (x >= left && x < right && y >= top && y < bottom);
}

bool isInside(int x, int y, int left, int top, int right, int bottom)
{
    return (x >= left && x < right && y >= top && y < bottom);
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
    Result result{0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}};

    // 1. Pre-compute

    // 1.1. Define pre-computed vectors
    vector<double> avgColor = vector<double>(N_COLOR_CHAN * nx * ny, 0.);
    avgColor[0] = data[0];
    avgColor[1] = data[1];
    avgColor[2] = data[2];

    // 1.1. Pre-compute 3 color channels
    // NOTE: HoangLe [May-21]: Can use multithread, z-order here
    for (int xy = 1; xy < ny * nx; ++xy) // Start at 2 because position [0, 0] is pre-assigned above
    {
        int x, y;
        xy2each(xy, nx, x, y);

        for (int c = 0; c < N_COLOR_CHAN; ++c)
        {
            double X = avgColor[c + 3 * x + 3 * nx * max(y - 1, 0)];
            double Y = avgColor[c + 3 * max(x - 1, 0) + 3 * nx * y];
            double Z = avgColor[c + 3 * max(x - 1, 0) + 3 * nx * max(y - 1, 0)];
            double current = data[c + 3 * x + 3 * nx * y];
            int Xn = getNumCells(x, y - 1, false);
            int Yn = getNumCells(x - 1, y, false);
            int Zn = getNumCells(x - 1, y - 1, false);
            int n = getNumCells(x, y, false);

            avgColor[c + 3 * x + 3 * nx * y] = 1.0 / n * (X * Xn + Y * Yn - Z * Zn + current);
        }
    }

    // 1.2. Pre-compute error
    // 1.2.1. Calculate possible pairs of topleft-bottomright
    vector<TLBR> pairs_tlbr = vector<TLBR>(ny * nx * ny * nx);
    int nPairsValid = 0;
    // NOTE: HoangLe [May-21]: Can use multithread
    for (int br = nx + 1; br <= (ny + 1) * (nx + 1); ++br)
    {
        int right, bottom;
        xy2each(br, nx + 1, right, bottom);

        for (int top = 0; top < bottom; ++top)
            for (int left = 0; left < right; ++left)
            {
                pairs_tlbr[nPairsValid] = make_tuple(top, left, bottom, right);
                ++nPairsValid;
            }
    }

    // 1.2.2. Pre-compute error with given valid tl-br pairs
    vector<double> errors = vector<double>(nPairsValid, 0);
    vector<double> inner = vector<double>(N_COLOR_CHAN * nPairsValid, 0);
    vector<double> outer = vector<double>(N_COLOR_CHAN * nPairsValid, 0);

// NOTE: HoangLe [May-21]: Can use multithread here
#pragma omp parallel for schedule(static, 1)
    for (int tlbr = 0; tlbr < nPairsValid; ++tlbr)
    {
        int left, top, right, bottom;
        top = get<0>(pairs_tlbr[tlbr]);
        left = get<1>(pairs_tlbr[tlbr]);
        bottom = get<2>(pairs_tlbr[tlbr]);
        right = get<3>(pairs_tlbr[tlbr]);

        // Determine outer[c] and inner[c]
        for (int c = 0; c < N_COLOR_CHAN; ++c)
        {
            for (int y = 0; y < ny; ++y)
                for (int x = 0; x < nx; ++x)
                {
                    if (isInside(x, y, left, top, right, bottom))
                    {
                        inner[tlbr * N_COLOR_CHAN + c] += data[c + 3 * x + 3 * nx * y];
                    }
                    else
                    {
                        outer[tlbr * N_COLOR_CHAN + c] += data[c + 3 * x + 3 * nx * y];
                    }
                }

            inner[tlbr * N_COLOR_CHAN + c] /= getNumCells(left, top, right, bottom);
            outer[tlbr * N_COLOR_CHAN + c] /= getNumCellsOutside(left, top, right, bottom, nx, ny);
        }

        // Determine the error
        for (int y = 0; y < ny; ++y)
            for (int x = 0; x < nx; ++x)
                for (int c = 0; c < N_COLOR_CHAN; ++c)
                {
                    if (isInside(x, y, left, top, right, bottom))
                    {
                        errors[tlbr] += pow(inner[tlbr * N_COLOR_CHAN + c] - data[c + 3 * x + 3 * nx * y], 2);
                    }
                    else
                    {
                        errors[tlbr] += pow(outer[tlbr * N_COLOR_CHAN + c] - data[c + 3 * x + 3 * nx * y], 2);
                    }
                }
    }

    // 2. Choose best tl, br
    // NOTE: HoangLe [May-21]: Can use multithread here
    double bestError = inf;
    for (int tlbr = 0; tlbr < nPairsValid; ++tlbr)
    {
        int left, top, right, bottom;
        top = get<0>(pairs_tlbr[tlbr]);
        left = get<1>(pairs_tlbr[tlbr]);
        bottom = get<2>(pairs_tlbr[tlbr]);
        right = get<3>(pairs_tlbr[tlbr]);

        if (bestError > errors[tlbr])
        {
            bestError = errors[tlbr];

            result.x0 = left;
            result.y0 = top;
            result.x1 = right;
            result.y1 = bottom;
            result.outer[0] = static_cast<float>(outer[N_COLOR_CHAN * tlbr + 0]);
            result.outer[1] = static_cast<float>(outer[N_COLOR_CHAN * tlbr + 1]);
            result.outer[2] = static_cast<float>(outer[N_COLOR_CHAN * tlbr + 2]);
            result.inner[0] = static_cast<float>(inner[N_COLOR_CHAN * tlbr + 0]);
            result.inner[1] = static_cast<float>(inner[N_COLOR_CHAN * tlbr + 1]);
            result.inner[2] = static_cast<float>(inner[N_COLOR_CHAN * tlbr + 2]);
        }
    }

    // Testing
    // Testing::printColor<vector<double>>(avgColor, nx, ny);

    return result;
}

#include <cmath>
#include <vector>

using namespace std;

constexpr int N_PER_PACK = 4;

typedef double double4 __attribute__((vector_size(N_PER_PACK * sizeof(double))));

const double4 vZeros = {0., 0., 0., 0.};

static inline double sumInternal(double4 x)
{
    double sum = 0;
    for (int i = 0; i < N_PER_PACK; ++i)
        sum += x[i];

    return sum;
}

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/
void correlate(int ny, int nx, const float *data, float *result)
{
    vector<double> norm(ny * nx, 0.);

    // 1. row-wise 0-mean normalization

#pragma omp parallel for
    for (int i = 0; i < ny; ++i)
    {
        double mean = 0;
        for (int j = 0; j < nx; ++j)
            mean += data[i * nx + j];

        mean /= nx;

        for (int j = 0; j < nx; ++j)
            norm[i * nx + j] = data[i * nx + j] - mean;
    }

    // 2. row-wise square-sum normalization

#pragma omp parallel for
    for (int i = 0; i < ny; ++i)
    {
        double sq_sum = 0;
        for (int j = 0; j < nx; ++j)
        {
            sq_sum += norm[i * nx + j] * norm[i * nx + j];
        }
        sq_sum = sqrt(sq_sum);
        for (int j = 0; j < nx; ++j)
            norm[i * nx + j] = norm[i * nx + j] / sq_sum;
    }

    // 3. upper-triangular matmul
    // 3.1. Expand vector 'norm'
    constexpr int bs = 3;
    int nPacks = (nx + N_PER_PACK - 1) / N_PER_PACK;
    int nyPadded = (ny + bs - 1) / bs * bs;

    // 3.2. Move 'norm' to new padded
    vector<double4> vnorm(nyPadded * nPacks, vZeros);
#pragma omp parallel for
    for (int i = 0; i < ny; ++i)
        for (int j = 0; j < nx; ++j)
            vnorm[nPacks * i + j / N_PER_PACK][j % N_PER_PACK] = norm[nx * i + j];

// 3.3. Start calculating
#pragma omp parallel for
    for (int i = 0; i < ny; ++i)
        result[ny * i + i] = 1.0;

#pragma omp parallel for
    for (int i1 = 0; i1 < nyPadded - bs; i1 += bs)
        for (int i2 = i1 + bs; i2 < nyPadded; i2 += bs)
        {
            double4 vtmp[bs][bs];

            // Initialize
            for (int i = 0; i < bs; ++i)
                for (int j = 0; j < bs; ++j)
                    vtmp[i][j] = vZeros;

            // Looping
            for (int k = 0; k < nPacks; ++k)
            {
                double4 v1 = vnorm[nPacks * (i1 + 0) + k];
                double4 v2 = vnorm[nPacks * (i1 + 1) + k];
                double4 v3 = vnorm[nPacks * (i1 + 2) + k];
                double4 v4 = vnorm[nPacks * (i2 + 0) + k];
                double4 v5 = vnorm[nPacks * (i2 + 1) + k];
                double4 v6 = vnorm[nPacks * (i2 + 2) + k];

                vtmp[0][0] = v1 * v4 + vtmp[0][0];
                vtmp[0][1] = v1 * v5 + vtmp[0][1];
                vtmp[0][2] = v1 * v6 + vtmp[0][2];
                vtmp[1][0] = v2 * v4 + vtmp[1][0];
                vtmp[1][1] = v2 * v5 + vtmp[1][1];
                vtmp[1][2] = v2 * v6 + vtmp[1][2];
                vtmp[2][0] = v3 * v4 + vtmp[2][0];
                vtmp[2][1] = v3 * v5 + vtmp[2][1];
                vtmp[2][2] = v3 * v6 + vtmp[2][2];
            }

            // Assign back to 'result'
            for (int i = 0; i < bs; ++i)
                for (int j = 0; j < bs; ++j)
                    if (i1 + i < ny && i2 + j < ny)
                        result[ny * (i1 + i) + i2 + j] = (float)sumInternal(vtmp[i][j]);
        }

#pragma omp parallel for
    for (int i1 = 0; i1 < nyPadded; i1 += bs)
    {
        double4 vtmp[] = {vZeros, vZeros, vZeros};
        for (int k = 0; k < nPacks; ++k)
        {
            double4 v1 = vnorm[nPacks * (i1 + 0) + k];
            double4 v2 = vnorm[nPacks * (i1 + 1) + k];
            double4 v3 = vnorm[nPacks * (i1 + 2) + k];

            vtmp[0] += v1 * v2;
            vtmp[1] += v1 * v3;
            vtmp[2] += v2 * v3;
        }
        if (i1 + 1 < ny)
            result[ny * (i1 + 0) + (i1 + 1)] = (float)sumInternal(vtmp[0]);
        if (i1 + 2 < ny)
        {
            result[ny * (i1 + 0) + (i1 + 2)] = (float)sumInternal(vtmp[1]);
            result[ny * (i1 + 1) + (i1 + 2)] = (float)sumInternal(vtmp[2]);
        }
    }
}
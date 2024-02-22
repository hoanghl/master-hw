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
    constexpr int bs = 8;
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
    for (int i1 = 0; i1 < nyPadded; i1 += bs)
        for (int i2 = i1; i2 < nyPadded; i2 += bs)
        {
            double4 vtmp[bs][bs];

            // Initialize
            for (int i = 0; i < bs; ++i)
                for (int j = 0; j < bs; ++j)
                    vtmp[i][j] = vZeros;

            // Looping
            for (int k = 0; k < nPacks; ++k)
            {
                double4 v1_1 = vnorm[nPacks * (i1 + 0) + k];
                double4 v1_2 = vnorm[nPacks * (i1 + 1) + k];
                double4 v1_3 = vnorm[nPacks * (i1 + 2) + k];
                double4 v1_4 = vnorm[nPacks * (i1 + 3) + k];
                double4 v1_5 = vnorm[nPacks * (i1 + 4) + k];
                double4 v1_6 = vnorm[nPacks * (i1 + 5) + k];
                double4 v1_7 = vnorm[nPacks * (i1 + 6) + k];
                double4 v1_8 = vnorm[nPacks * (i1 + 7) + k];

                double4 v2_1 = vnorm[nPacks * (i2 + 0) + k];
                double4 v2_2 = vnorm[nPacks * (i2 + 1) + k];
                double4 v2_3 = vnorm[nPacks * (i2 + 2) + k];
                double4 v2_4 = vnorm[nPacks * (i2 + 3) + k];
                double4 v2_5 = vnorm[nPacks * (i2 + 4) + k];
                double4 v2_6 = vnorm[nPacks * (i2 + 5) + k];
                double4 v2_7 = vnorm[nPacks * (i2 + 6) + k];
                double4 v2_8 = vnorm[nPacks * (i2 + 7) + k];

                vtmp[0][0] = v1_1 * v2_1 + vtmp[0][0];
                vtmp[0][1] = v1_1 * v2_2 + vtmp[0][1];
                vtmp[0][2] = v1_1 * v2_3 + vtmp[0][2];
                vtmp[0][3] = v1_1 * v2_4 + vtmp[0][3];
                vtmp[0][4] = v1_1 * v2_5 + vtmp[0][4];
                vtmp[0][5] = v1_1 * v2_6 + vtmp[0][5];
                vtmp[0][6] = v1_1 * v2_7 + vtmp[0][6];
                vtmp[0][7] = v1_1 * v2_8 + vtmp[0][7];

                vtmp[1][0] = v1_2 * v2_1 + vtmp[1][0];
                vtmp[1][1] = v1_2 * v2_2 + vtmp[1][1];
                vtmp[1][2] = v1_2 * v2_3 + vtmp[1][2];
                vtmp[1][3] = v1_2 * v2_4 + vtmp[1][3];
                vtmp[1][4] = v1_2 * v2_5 + vtmp[1][4];
                vtmp[1][5] = v1_2 * v2_6 + vtmp[1][5];
                vtmp[1][6] = v1_2 * v2_7 + vtmp[1][6];
                vtmp[1][7] = v1_2 * v2_8 + vtmp[1][7];

                vtmp[2][0] = v1_3 * v2_1 + vtmp[2][0];
                vtmp[2][1] = v1_3 * v2_2 + vtmp[2][1];
                vtmp[2][2] = v1_3 * v2_3 + vtmp[2][2];
                vtmp[2][3] = v1_3 * v2_4 + vtmp[2][3];
                vtmp[2][4] = v1_3 * v2_5 + vtmp[2][4];
                vtmp[2][5] = v1_3 * v2_6 + vtmp[2][5];
                vtmp[2][6] = v1_3 * v2_7 + vtmp[2][6];
                vtmp[2][7] = v1_3 * v2_8 + vtmp[2][7];

                vtmp[3][0] = v1_4 * v2_1 + vtmp[3][0];
                vtmp[3][1] = v1_4 * v2_2 + vtmp[3][1];
                vtmp[3][2] = v1_4 * v2_3 + vtmp[3][2];
                vtmp[3][3] = v1_4 * v2_4 + vtmp[3][3];
                vtmp[3][4] = v1_4 * v2_5 + vtmp[3][4];
                vtmp[3][5] = v1_4 * v2_6 + vtmp[3][5];
                vtmp[3][6] = v1_4 * v2_7 + vtmp[3][6];
                vtmp[3][7] = v1_4 * v2_8 + vtmp[3][7];

                vtmp[4][0] = v1_5 * v2_1 + vtmp[4][0];
                vtmp[4][1] = v1_5 * v2_2 + vtmp[4][1];
                vtmp[4][2] = v1_5 * v2_3 + vtmp[4][2];
                vtmp[4][3] = v1_5 * v2_4 + vtmp[4][3];
                vtmp[4][4] = v1_5 * v2_5 + vtmp[4][4];
                vtmp[4][5] = v1_5 * v2_6 + vtmp[4][5];
                vtmp[4][6] = v1_5 * v2_7 + vtmp[4][6];
                vtmp[4][7] = v1_5 * v2_8 + vtmp[4][7];

                vtmp[5][0] = v1_6 * v2_1 + vtmp[5][0];
                vtmp[5][1] = v1_6 * v2_2 + vtmp[5][1];
                vtmp[5][2] = v1_6 * v2_3 + vtmp[5][2];
                vtmp[5][3] = v1_6 * v2_4 + vtmp[5][3];
                vtmp[5][4] = v1_6 * v2_5 + vtmp[5][4];
                vtmp[5][5] = v1_6 * v2_6 + vtmp[5][5];
                vtmp[5][6] = v1_6 * v2_7 + vtmp[5][6];
                vtmp[5][7] = v1_6 * v2_8 + vtmp[5][7];

                vtmp[6][0] = v1_7 * v2_1 + vtmp[6][0];
                vtmp[6][1] = v1_7 * v2_2 + vtmp[6][1];
                vtmp[6][2] = v1_7 * v2_3 + vtmp[6][2];
                vtmp[6][3] = v1_7 * v2_4 + vtmp[6][3];
                vtmp[6][4] = v1_7 * v2_5 + vtmp[6][4];
                vtmp[6][5] = v1_7 * v2_6 + vtmp[6][5];
                vtmp[6][6] = v1_7 * v2_7 + vtmp[6][6];
                vtmp[6][7] = v1_7 * v2_8 + vtmp[6][7];

                vtmp[7][0] = v1_8 * v2_1 + vtmp[7][0];
                vtmp[7][1] = v1_8 * v2_2 + vtmp[7][1];
                vtmp[7][2] = v1_8 * v2_3 + vtmp[7][2];
                vtmp[7][3] = v1_8 * v2_4 + vtmp[7][3];
                vtmp[7][4] = v1_8 * v2_5 + vtmp[7][4];
                vtmp[7][5] = v1_8 * v2_6 + vtmp[7][5];
                vtmp[7][6] = v1_8 * v2_7 + vtmp[7][6];
                vtmp[7][7] = v1_8 * v2_8 + vtmp[7][7];
            }

            // Assign back to 'result'
            for (int i = 0; i < bs; ++i)
                for (int j = 0; j < bs; ++j)
                    if (i1 + i < ny && i2 + j < ny && i1 + i < i2 + j)
                        result[ny * (i1 + i) + i2 + j] = (float)sumInternal(vtmp[i][j]);
        }
}
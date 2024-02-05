// #include "cp.hpp"

#include <cmath>
#include <vector>

using namespace std;

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

    for (int i = 0; i < ny; ++i)
    {
        double sq_sum = 0;
        for (int j = 0; j < nx; ++j)
            sq_sum += pow(norm[i * nx + j], 2);
        sq_sum = sqrt(sq_sum);
        for (int j = 0; j < nx; ++j)
            norm[i * nx + j] = norm[i * nx + j] / sq_sum;
    }

    const int nCon = 4;
    int nxPadded = nx;
    if (nx % nCon != 0)
        nxPadded = (nx / nCon + 1) * nCon;

    // 3. Move to new array
    vector<double> normPadded(ny * nxPadded, 0.);
    for (int i = 0; i < ny; ++i)
        for (int j = 0; j < nx; ++j)
            normPadded[nxPadded * i + j] = norm[nx * i + j];

    // 3. upper-triangular matmul

    for (int i = 0; i < ny; ++i)
        for (int j = i; j < ny; ++j)
        {
            double corArr[nCon];
            for (int l = 0; l < nCon; ++l)
                corArr[l] = 0;
            for (int k = 0; k < nxPadded / nCon; ++k)
            {
                for (int l = 0; l < nCon; ++l)
                    corArr[l] += normPadded[i * nxPadded + k * nCon + l] * normPadded[j * nxPadded + k * nCon + l];
            }

            double cor = 0;
            for (int l = 0; l < nCon; ++l)
                cor += corArr[l];

            result[i * ny + j] = (float)cor;
        }
}

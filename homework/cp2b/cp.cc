#include <cmath>

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
    double *norm = new double[ny * nx];

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
            sq_sum += norm[i * nx + j] * norm[i * nx + j];
        sq_sum = sqrt(sq_sum);
        for (int j = 0; j < nx; ++j)
            norm[i * nx + j] = norm[i * nx + j] / sq_sum;
    }

    // 3. upper-triangular matmul

#pragma omp parallel for
    for (int i = 0; i < ny; ++i)
        for (int j = i; j < ny; ++j)
        {
            double cor = 0;
            for (int k = 0; k < nx; ++k)
                cor += norm[i * nx + k] * norm[j * nx + k];

            result[i * ny + j] = (float)cor;
        }
}

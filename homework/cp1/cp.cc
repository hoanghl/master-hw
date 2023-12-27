// #include "cp.hpp"

#include <math.h>
#include <iostream>
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
    vector<double> tmp(ny * nx, 0.);

    // 1. row-wise 0-mean normalization
    // cout << "Step 1" << endl;

    for (int i = 0; i < ny; ++i)
    {
        double mean = 0;
        for (int j = 0; j < nx; ++j)
            mean += data[i * nx + j];

        mean /= nx;

        for (int j = 0; j < nx; ++j)
            tmp[i * nx + j] = data[i * nx + j] - mean;
    }

    // 2. row-wise square-sum normalization
    // cout << "Step 2" << endl;

    for (int i = 0; i < ny; ++i)
    {
        double sq_sum = 0;
        for (int j = 0; j < nx; ++j)
        {
            sq_sum += tmp[i * nx + j] * tmp[i * nx + j];
        }
        sq_sum = sqrt(sq_sum);
        for (int j = 0; j < nx; ++j)
            tmp[i * nx + j] = tmp[i * nx + j] / sq_sum;
    }

    // 3. upper-triangular matmul
    // cout << "Step 3" << endl;

    for (int i = 0; i < ny; ++i)
        for (int j = i; j < ny; ++j)
        {
            double cor = 0;
            for (int k = 0; k < nx; ++k)
                cor += tmp[i * nx + k] * tmp[j * nx + k];

            result[i * ny + j] = (float)cor;
        }
}

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

    // 1. Pad nx so that nx is divisible by 4
    int nx_final = nx;
    if (nx % 4 != 0)
    {
        nx_final = (nx / 4 + 1) * 4;
    }
    int nx_para = nx_final / 4;

    vector<double> data_padded(ny * nx_final, 0.);
    vector<double> tmp(ny * nx_final, 0.);

    // 0. Move data from old to new
    for (int i = 0; i < ny; ++i)
        for (int j = 0; j < nx; ++j)
            data_padded[i * nx_final + j] = data[i * nx + j];

    // 1. row-wise 0-mean normalization
    // #pragma omp parallel for
    vector<double> means(4, 0.);
    for (int i = 0; i < ny; ++i)
    {
        // Calculate mean by folding technique
        for (int n_fold = 0; n_fold < 4; ++n_fold)
            for (int j = 0; j < nx_para; ++j)
                means[n_fold] += data_padded[i * nx_final + j * 4 + n_fold];

        double mean = 0;
        for (int n_fold = 0; n_fold < 4; ++n_fold)
        {
            mean += means[n_fold];
            means[n_fold] = 0;
        }
        mean /= nx;

        // Apply normalization
        for (int j = 0; j < nx; ++j)
            tmp[i * nx_final + j] = data_padded[i * nx_final + j] - mean;
    }

    // 2. row-wise square-sum normalization
    // cout << "Step 2" << endl;

    // #pragma omp parallel for
    vector<double> sq_sums(4, 0.);
    for (int i = 0; i < ny; ++i)
    {

        // Calculate mean by folding technique
        for (int n_fold = 0; n_fold < 4; ++n_fold)
            for (int j = 0; j < nx_para; ++j)
                sq_sums[n_fold] += tmp[i * nx_final + j * 4 + n_fold] * tmp[i * nx_final + j * 4 + n_fold];
        double sq_sum = 0;
        for (int n_fold = 0; n_fold < 4; ++n_fold)
        {
            sq_sum += sq_sums[n_fold];
            sq_sums[n_fold] = 0;
        }
        sq_sum = sqrt(sq_sum);

        for (int j = 0; j < nx; ++j)
            tmp[i * nx_final + j] = tmp[i * nx_final + j] / sq_sum;
    }

    // 3. upper-triangular matmul
    // cout << "Step 3" << endl;

    // #pragma omp parallel for
    vector<double> cors(4, 0.);
    for (int i = 0; i < ny; ++i)
        for (int j = i; j < ny; ++j)
        {

            // Calculate mean by folding technique
            for (int n_fold = 0; n_fold < 4; ++n_fold)
                for (int k = 0; k < nx_para; ++k)
                    cors[n_fold] += tmp[i * nx_final + k * 4 + n_fold] * tmp[j * nx_final + k * 4 + n_fold];
            double cor = 0;
            for (int n_fold = 0; n_fold < 4; ++n_fold)
            {
                cor += cors[n_fold];
                cors[n_fold] = 0;
            }

            result[i * ny + j] = (float)cor;
        }
}
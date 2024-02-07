#include <cmath>
#include <iomanip>
#include <iostream>
#include <time.h>

#include "ex4.hpp"

using namespace std;

const int N = 1000;

int main(int argc, char const *argv[])
{
    // Initializa
    double **A = genMatrix<double>(N);
    double **B = genMatrix<double>(N);
    double **C = genMatrixEmpty<double>(N);

    clock_t start = clock();
#pragma omp parallel for
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k)
                C[i][j] += A[k][i] * B[k][j];
    clock_t end = clock();
    double elapsed = double(end - start) / CLOCKS_PER_SEC;
    cout << "CPU time (second): " << setprecision(12) << elapsed << endl;
    return 0;
}

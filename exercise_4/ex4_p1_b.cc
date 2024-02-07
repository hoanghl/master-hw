#include <iostream>
#include <random>
#include <iomanip>
#include <time.h>

#include "ex4.hpp"

using namespace std;

void probB_Original(int n, float **a, float **b, float *c)
{
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            a[i][j] = b[i][j] / c[i];
}

void probB_Optimized(int n, float **a, float **b, float *c)
{
#pragma omp parallel for
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            a[i][j] = b[i][j] / c[i];
}

int main(int argc, char const *argv[])
{
    cout << "* Problem A: " << endl;

    int n = 3000;
    float **a = new float *[n], **b = new float *[n], *c = new float[n];

    // Fill a, b and c
    for (int i = 0; i < n; ++i)
    {
        a[i] = new float[n];
        b[i] = new float[n];

        for (int j = 0; j < n; ++j)
        {
            a[i][j] = getRand<float>();
            b[i][j] = getRand<float>();
        }

        c[i] = getRand<float>();
    }

    cout << "= Original ==============" << endl;

    clock_t start = clock();
    probB_Original(n, a, b, c);
    clock_t end = clock();
    double elapsed = double(end - start) / CLOCKS_PER_SEC;
    cout << "CPU time (second): " << setprecision(12) << elapsed << endl;

    cout << "= Optimized =============" << endl;

    start = clock();
    probB_Optimized(n, a, b, c);
    end = clock();
    elapsed = double(end - start) / CLOCKS_PER_SEC;
    cout << "CPU time (second): " << setprecision(12) << elapsed << endl;

    return 0;
}

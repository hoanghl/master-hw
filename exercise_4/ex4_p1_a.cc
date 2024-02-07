#include <iostream>
#include <random>
#include <iomanip>
#include <time.h>

#include "ex4.hpp"

using namespace std;

void probA_Original(int n, float *a, float *b)
{
    for (int i = 0; i < n - 1; ++i)
    {
        if (i < 500)
            a[i] = 4.0 * b[i] + b[i + 1];
        else
            a[i] = 4.0 * b[i + 1] + b[i];
    }
}

void probA_Optimized(int n, float *a, float *b)
{
#pragma omp parallel for
    for (int i = 0; i < 500; ++i)
        a[i] = 4.0 * b[i] + b[i + 1];

#pragma omp parallel for
    for (int i = 500; i <= n - 1; ++i)
        a[i] = 4.0 * b[i + 1] + b[i];
}

int main(int argc, char const *argv[])
{
    cout << "* Problem A: " << endl;

    int n = 10000000;
    float *a = new float[n], *b = new float[n];

    // Fill b
    for (int i = 0; i < n; ++i)
    {
        b[i] = getRand<float>();
    }

    cout << "= Original ==============" << endl;

    clock_t start = clock();
    probA_Original(n, a, b);
    clock_t end = clock();
    double elapsed = double(end - start) / CLOCKS_PER_SEC;
    cout << "CPU time (second): " << setprecision(12) << elapsed << endl;

    cout << "= Optimized =============" << endl;

    start = clock();
    probA_Optimized(n, a, b);
    end = clock();
    elapsed = double(end - start) / CLOCKS_PER_SEC;
    cout << "CPU time (second): " << setprecision(12) << elapsed << endl;

    return 0;
}

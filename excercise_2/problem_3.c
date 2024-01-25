#include <stdio.h>
#include <math.h>
#include <time.h>
#define DENO 1000000

int main()
{
    unsigned int n = 100000000;

    clock_t t1, t2;
    double cputime;

    // Calculate with: float
    t1 = clock();

    float result = 0;
    for (unsigned int k = 0; k <= n; k++)
        result += expf(sinf(k / DENO));

    t2 = clock();
    cputime = (double)(t2 - t1) / CLOCKS_PER_SEC;
    printf("= Case: Float ==========================\n");
    printf("- Result: %f\n", result);
    printf("- CPU time in seconds: %g\n", cputime);
    printf("========================================\n");

    // Calculate with: double
    t1 = clock();

    double result_d = 0;
    for (unsigned int k = 0; k <= n; k++)
        result_d += exp(sin(k / DENO));

    t2 = clock();
    cputime = (double)(t2 - t1) / CLOCKS_PER_SEC;
    printf("= Case: Double =========================\n");
    printf("- Result: %f\n", result_d);
    printf("- CPU time in seconds: %g\n", cputime);
    printf("========================================\n");

    // Calculate with: 128
    t1 = clock();

    __float128 result_128 = 0;
    for (unsigned int k = 0; k <= n; k++)
        result_128 += expl(sinl(k / DENO));

    t2 = clock();
    result_d = (double)result_128;
    cputime = (double)(t2 - t1) / CLOCKS_PER_SEC;
    printf("= Case: Double =========================\n");
    printf("- Result: %f\n", result_d);
    printf("- CPU time in seconds: %g\n", cputime);
    printf("========================================\n");
}
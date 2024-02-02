#include <iostream>
#include "math.h"
#include <time.h>
#include <sys/time.h>
#include <fstream>
#include <iomanip>

#define N 50000

using namespace std;

double getRand()
{
    return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}

void fn_measure()
{
    double arr[N];
    arr[0] = getRand();

    for (int i = 1; i < N; ++i)
    {
        for (int j = i - 1; j >= 0; --j)
        {
            double a = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
            double b = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);

            arr[j] = a * b - sqrt(abs(arr[i])) * sqrt(abs(a - b));
        }
    }

    for (int k = 0; k < 30; ++k)
    {
        for (int i = 0; i < N; ++i)
            cout << setprecision(10) << arr[i] << endl;
    }
}

int main()
{
    // Measure
    clock_t start = clock();
    fn_measure();
    clock_t end = clock();
    double elapsed = double(end - start) / CLOCKS_PER_SEC;
    cout << "CPU time: " << setprecision(5) << elapsed << endl;

    struct timeval beginWall, endWall;
    gettimeofday(&beginWall, 0);
    fn_measure();
    gettimeofday(&endWall, 0);
    long seconds = endWall.tv_sec - beginWall.tv_sec;
    long microseconds = endWall.tv_usec - beginWall.tv_usec;
    elapsed = seconds + microseconds * 1e-6;
    cout << "Wall time: " << setprecision(5) << elapsed;
}
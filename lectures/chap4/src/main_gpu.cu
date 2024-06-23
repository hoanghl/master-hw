#include <chrono>
#include <cstdlib>
#include <iostream>

#include "gpu.hpp"
#include "utils.hpp"

using namespace std;
using namespace std::chrono;

int main(int argc, char const *argv[])
{
    int n = 6000;

    // 1. Create data
    printf("1. Create data\n");

    float *data = Data::genData(n);
    // Input *d = Data::readData(Data::FILE1);
    // if (d == nullptr)
    //     return -1;
    // n = d->ny;
    // float *data = d->d;

    // 2. Start
    printf("2. Start with ny = %d, nx = %d\n", n, n);

    float *result = new float[n * n];

    printf("= c4v0\n");

    auto t1 = high_resolution_clock::now();
    GPU::c4v0(n, data, result);
    auto t2 = high_resolution_clock::now();

    duration<double, std::milli> ms_double = t2 - t1;
    printf("=> Running time: %12.6g\n", ms_double.count());

    printf("= c4v1\n");

    t1 = high_resolution_clock::now();
    GPU::c4v1(n, data, result);
    t2 = high_resolution_clock::now();

    ms_double = t2 - t1;
    printf("=> Running time: %12.6g\n", ms_double.count());

    // delete[] d->d;
    delete data;
    delete[] result;

    return 0;
}
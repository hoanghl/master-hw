#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <vector>

#include "exp.hpp"

using namespace std;
using namespace std::chrono;

typedef unsigned long long data_t;
// typedef int data_t;
constexpr int NUM_BASE = 50;

void test(int n, data_t *&data)
{
    for (int i = 0; i < n; ++i)
    {
        printf("i = %3d: %3llu\n", i, data[i]);
    }
}

void test(int n, vector<data_t> &data)
{
    for (int i = 0; i < n; ++i)
    {
        printf("i = %3d: %3llu\n", i, data[i]);
    }
}

void test(int l, int r, data_t *&data)
{
    for (int i = l; i < r; ++i)
    {
        printf("i = %3d: %3llu\n", i, data[i]);
    }
}

void mergesort_merge(data_t *data, vector<data_t> &aux, int l, int m, int r)
{
    int i1 = l, i2 = m;
    for (int i = l; i < r; ++i)
    {
        if (i2 >= r)
        {
            aux[i] = data[i1];
            ++i1;
        }
        else if (i1 >= m)
        {
            aux[i] = data[i2];
            ++i2;
        }
        else if (data[i1] < data[i2])
        {
            aux[i] = data[i1];
            ++i1;
        }
        else
        {
            aux[i] = data[i2];
            ++i2;
        }
    }

    for (int i = l; i < r; ++i)
    {
        data[i] = aux[i];
    }
}

void mergesort_sort(data_t *data, vector<data_t> &aux, int l, int r)
{
    if (r - l <= NUM_BASE)
    {
        sort(data + l, data + r);
        return;
    }

    int m = (r - l) / 2 + l;

    mergesort_sort(data, aux, l, m);
    mergesort_sort(data, aux, m, r);
    mergesort_merge(data, aux, l, m, r);
}

void mergesort_sort2(data_t *data, vector<data_t> &aux, int l, int r)
{
    if (r - l <= NUM_BASE)
    {
        sort(data + l, data + r);
        return;
    }

    int m = (r - l) / 2 + l;

#pragma omp task
    mergesort_sort(data, aux, l, m);
#pragma omp task
    mergesort_sort(data, aux, m, r);

    mergesort_merge(data, aux, l, m, r);
}

/////////////////////////////////////////////////////////////////////
void mergesort(int N, data_t *data)
{
    if (N <= NUM_BASE)
    {
        sort(data, data + N);
        return;
    }

    auto t1 = high_resolution_clock::now();

    vector<data_t> aux = vector<data_t>(N);
    int numEach = 0;

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int nT = omp_get_num_threads();
        numEach = (N + nT - 1) / nT;

        int left = numEach * tid, right = min(N, numEach * (tid + 1));

        if (left < right)
        {
            // printf("nT=%2d - numEach=%2d - tID=%2d - left=%2d - right=%2d\n", nT, numEach, tid, left, right);
            mergesort_sort(data, aux, left, right);
        }
    }

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    printf("1. Running time: %8.4f\n", ms_double.count());

    // printf("after each-seg-sort--\n");
    // test(N, data);
    // printf("after each-seg-sort-end--\n");

    // Merge all segments
    t1 = high_resolution_clock::now();

    int lastNumEach = numEach;
    numEach = min(numEach * 2, N);

    while (numEach <= N)
    {
#pragma omp parallel for
        for (int i = 0; i < N; i += numEach)
        {

            int i1 = i, i2 = i + lastNumEach;
            int lim_i1 = min(N, i1 + lastNumEach), lim_i2 = min(N, i + lastNumEach * 2);

            for (int j = i; j < min(N, i + numEach); ++j)
            {
                if (i2 >= lim_i2)
                {
                    aux[j] = data[i1];
                    ++i1;
                }
                else if (i1 >= lim_i1)
                {
                    aux[j] = data[i2];
                    ++i2;
                }
                else if (data[i1] < data[i2])
                {
                    aux[j] = data[i1];
                    ++i1;
                }
                else
                {
                    aux[j] = data[i2];
                    ++i2;
                }
            }

            // for (int j = i; j < min(N, i + numEach); ++j)
            //     data[j] = aux[j];
        }

// Update back to data
#pragma omp parallel for
        for (int i = 0; i < N; ++i)
            data[i] = aux[i];

        if (numEach == N)
            break;
        lastNumEach = numEach;
        numEach = min(numEach * 2, N);
    }

    t2 = high_resolution_clock::now();
    ms_double = t2 - t1;
    printf("2. Running time: %8.4f\n", ms_double.count());
}

void mergesort2(int N, data_t *data)
{
    if (N <= NUM_BASE)
    {
        sort(data, data + N);
        return;
    }

    int numEach = NUM_BASE;
    int lastNumEach = numEach;

    // Sort with base

    // auto t1 = high_resolution_clock::now();

#pragma omp parallel for
    for (int i = 0; i < N; i += numEach)
    {
        int interval = min(numEach, N - i);
        sort(data + i, data + i + interval);
    }

    // auto t2 = high_resolution_clock::now();
    // duration<double, std::milli> ms_double = t2 - t1;
    // printf("1. Running time: %8.4f\n", ms_double.count());

    // Start merging

    printf("1.\n");
    test(N, data);
    printf("1-End.\n");

    // t1 = high_resolution_clock::now();

    vector<data_t> aux = vector<data_t>(N);

    lastNumEach = numEach;
    numEach = min(numEach * 2, N);
    while (numEach <= N)
    {
        // #pragma omp parallel for
        for (int i = 0; i < N; i += numEach)
        {
            int i1 = i, i2 = i + lastNumEach;
            int lim_i1 = min(N, i1 + lastNumEach), lim_i2 = min(N, i + lastNumEach * 2);

            for (int j = i; j < min(N, i + numEach); ++j)
            {
                if (i2 >= lim_i2)
                {
                    aux[j] = data[i1];
                    ++i1;
                }
                else if (i1 >= lim_i1)
                {
                    aux[j] = data[i2];
                    ++i2;
                }
                else if (data[i1] < data[i2])
                {
                    aux[j] = data[i1];
                    ++i1;
                }
                else
                {
                    aux[j] = data[i2];
                    ++i2;
                }
            }

            // for (int j = i; j < min(N, i + numEach); ++j)
            //     data[j] = aux[j];
            int length = min(numEach, N - i);
            copy_n(aux.begin() + i, length, data + i);
        }

        printf("numEach=%d\n", numEach);
        test(N, data);
        printf("numEach=%d -End.\n", numEach);

        // Update back to data

        // #pragma omp parallel for
        // for (int i = 0; i < N; ++i)
        //     data[i] = aux[i];
        // copy(data, data + N, aux.begin());

        if (numEach == N)
            break;
        lastNumEach = numEach;
        numEach = min(numEach * 2, N);
    }

    // t2 = high_resolution_clock::now();
    // ms_double = t2 - t1;
    // printf("2. Running time: %8.4f\n", ms_double.count());
}

void mergesort3(int N, data_t *data)
{
    if (N <= NUM_BASE)
    {
        sort(data, data + N);
        return;
    }

    vector<data_t> aux = vector<data_t>(N);

#pragma omp parallel
    {
#pragma omp single
        mergesort_sort2(data, aux, 0, N);
    }
}

int main(int argc, char const *argv[])
{

    // Generate random vectors
    // constexpr int N = 101;
    // data_t *data = new data_t[N];
    // DataGen::genRand(N, data);

    Input *input = DataGen::readFile(FILE1);
    int N = input->n;
    data_t *data = input->data;

    printf("=> Before\n");
    test(N, data);

    // auto t1 = high_resolution_clock::now();
    mergesort2(N, data);
    // auto t2 = high_resolution_clock::now();
    // duration<double, std::milli> ms_double = t2 - t1;
    // printf("Running time: %8.4f\n", ms_double.count());

    printf("=> After\n");
    test(N, data);

    delete[] input->data;
    delete input;

    return 0;
}

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <vector>

using namespace std;

typedef unsigned long long data_t;
constexpr int NUM_BASE = 20;

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

void psort(int N, data_t *data)
{
    if (N <= NUM_BASE)
    {
        sort(data, data + N);
        return;
    }

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

    // Merge all segments

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

            // int length = min(numEach, N - i);
            // copy_n(aux.begin() + i, length, data + i);
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
}

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

using namespace std;

typedef unsigned long long data_t;

void psort(int N, data_t *data)
{
    // printf("N = %d\n", N);
    if (N < 5)
    {
        sort(data, data + N);
    }
    else
    {
        int numBase = 5;

        int numEach = numBase;
        int lastNumEach = numEach;

        // Sort with base
        for (int i = 0; i < N; i += numEach)
        {
            int interval = min(numEach, N - i);
            sort(data + i, data + i + interval);
        }

        // Start merging
        vector<data_t> aux = vector<data_t>(N);

        lastNumEach = numEach;
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
}

#include <random>
#include <iostream>
#include <time.h>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <chrono>

using namespace std;
using namespace std::chrono;

const uint N = 10000;

int **readFile(string fileName)
{
    int **d = new int *[N];
    for (uint i = 0; i < N; i++)
        d[i] = new int[N];

    fstream file(fileName, ios::in);
    if (file.is_open())
    {
        string s;
        uint i = 0, j = 0;
        while (getline(file, s))
        {
            d[i][j] = stoi(s);
            i += (j + 1) / N;
            j = (j + 1) % N;

            if (i == N)
                break;
        }

        file.close();
    }

    return d;
}

void writeFile(string fileName)
{
    srand(time(0));

    fstream file(fileName, ios::out);
    if (file.is_open())
    {
        for (uint i = 0, j = 0; i < N * N; ++i)
        {
            if (i == j * j)
            {
                ++j;
                file << 0 << '\n';
            }
            else
            {
                file << rand() % 300 << '\n';
            }
        }

        file.close();
    }
}

bool comp(int a, int b)
{
    if (a <= b)
        return a;
    else
        return b;
}
void stepV1(int *r, const int *d, int n)
{
    auto start = high_resolution_clock::now();
    // #pragma omp parallel for
    for (uint i = 0; i < n; ++i)
        for (uint j = 0; j < n; ++j)
        {
            int v = 100000;
            for (uint k = 0; k < n; ++k)
            {
                int x = d[n * i + k];
                int y = d[n * k + j];
                int z = x + y;
                v = std::min(v, z);
            }
            r[n * i + j] = v;
        }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);

    cout << "Execution time: " << duration.count();
}

void stepV2(int *r, const int *d, int n)
{
    std::vector<int> t(n * n);
#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            t[n * j + i] = d[n * i + j];
        }
    }

    auto start = high_resolution_clock::now();

#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            int v = 100000;
            for (int k = 0; k < n; ++k)
            {
                int x = d[n * i + k];
                int y = t[n * j + k];
                int z = x + y;
                v = std::min(v, z);
            }
            r[n * i + j] = v;
        }
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);

    cout << "Execution time: " << duration.count();
}

int main()
{
    // writeFile("test1.txt");
    int **d = readFile("test1.txt");
    int *nums = new int[N * N];
    for (uint i = 0; i < N; i++)
        for (uint j = 0; j < N; j++)
            nums[i * N + j] = d[i][j];

    int *results = new int[N * N];

    stepV1(results, nums, N);
    // stepV2(results, nums, N);

    return 0;
}
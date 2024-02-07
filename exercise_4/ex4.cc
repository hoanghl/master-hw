#include "ex4.hpp"

#include <cmath>
#include <random>

using namespace std;

template <class T>
T getRand()
{
    return static_cast<T>(rand()) + static_cast<T>(rand()) / static_cast<T>(RAND_MAX);
}
template double getRand();
template float getRand();

template <class T>
T **genMatrix(int n)
{
    T **arr = new T *[n];
    for (int i = 0; i < n; ++i)
    {
        arr[i] = new T[n];
        for (int j = 0; j < n; ++j)
            arr[i][j] = getRand<T>();
    }

    return arr;
}

template double **genMatrix(int n);

template <class T>
T **genMatrixEmpty(int n)
{
    T **arr = new T *[n];
    for (int i = 0; i < n; ++i)
    {
        arr[i] = new T[n];
    }

    return arr;
}

template double **genMatrixEmpty(int n);

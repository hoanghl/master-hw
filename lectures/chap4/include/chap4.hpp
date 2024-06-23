#include <string>

#ifndef CHAP4_HPP
#define CHAP4_HPP

using namespace std;

typedef struct Input
{
    int ny, nx;
    float *d;
} Input;

namespace Data
{
    const string FILE1 = "data/1.txt";

    float *genData(int n);
    Input *readData(string filename);
};

namespace CPU
{
    void baseline(int n, const float *data, float *result);
    void c2v7(int n, const float *data, float *result);
}

namespace GPU
{
    static inline void check(cudaError_t err, const char *context)
    {
        if (err != cudaSuccess)
        {
            std::cerr << "CUDA error: " << context << ": "
                      << cudaGetErrorString(err) << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    static inline int divup(int a, int b)
    {
        return (a + b - 1) / b;
    }

#define CHECK(x) check(x, #x)

    void c4v0(int n, const float *data, float *result);
    // void c4v1(int n, const float *data, float *result);
    // void c4v2(int n, const float *data, float *result);
    // void c4v3(int n, const float *data, float *result);
}

#endif
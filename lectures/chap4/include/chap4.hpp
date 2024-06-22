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

#endif
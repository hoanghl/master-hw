#include <cstdlib>
#include <iostream>
#include <string>

using namespace std;

#ifndef UTILS_HPP
#define UTILS_HPP

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

#endif
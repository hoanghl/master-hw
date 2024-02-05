#include <string>

using namespace std;

#ifndef CC_HPP
#define CC_HPP

const string FILENAME = "test1.txt";

struct result
{
    int nx, ny;
    float *d;
};

void correlate(int ny, int nx, const float *data, float *result);

#endif
#include "chap2.hpp"
#include <iostream>
#include <chrono>

using namespace std;
using namespace std::chrono;

void step(float *r, const float *d, int n)
{
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        {
            float v = std::numeric_limits<float>::infinity();
            for (int k = 0; k < n; ++k)
            {
                float x = d[n * i + k];
                float y = d[n * k + j];
                float z = x + y;
                v = std::min(v, z);
            }

            r[n * i + j] = v;
        }
}

int main()
{

    chap2_result *result = readFile(STD_FILENAME);

    cout << "Finish reading data!" << endl;

    float *r = new float[result->n * result->n];

    cout << "Start stepping!" << endl;

    auto t1 = high_resolution_clock::now();
    step(r, result->d, result->n);
    auto t2 = high_resolution_clock::now();

    duration<double, std::milli> ms_double = t2 - t1;
    cout << "Running time: " << ms_double.count() << endl;

    return 0;
}
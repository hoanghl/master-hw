#include "chap4.hpp"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>

using namespace std;

float *Data::genData(int n)
{
    float *arr = new float[n * n];

    for (int i = 0; i < n * n; ++i)
        arr[i] = static_cast<float>(rand() * 1.0 / RAND_MAX + rand());

    return arr;
}

inline bool isFileExisted(const string &name)
{
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

Input *Data::readData(string filename)
{
    Input *d = nullptr;

    // Check file existence
    if (!isFileExisted(filename))
    {
        printf("Err: File not existed: %s\n", filename);
        return d;
    }

    fstream file(filename, ios::in);

    if (file.is_open())
    {
        d = new Input;
        string line;
        int i = -1;
        while (getline(file, line))
        {
            istringstream iss(line);

            if (i == -1)
            {
                ++i;

                iss >> d->ny;
                iss >> d->nx;

                d->d = new float[d->ny * d->nx];
            }
            else
            {
                for (int j = 0; j < d->nx; ++j)
                {
                    iss >> d->d[i];
                    ++i;
                }
            }
        }

        file.close();
    }

    return d;
}
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "cp.hpp"

using namespace std;

result *readFile()
{
    result *res = NULL;

    fstream file(FILENAME, ios::in);
    if (file.is_open())
    {
        string line;
        int nth_line = 0;
        res = new result;
        while (getline(file, line))
        {
            istringstream is(line);
            if (nth_line == 0)
            {
                is >> res->ny;
            }
            else if (nth_line == 1)
            {
                is >> res->nx;
            }
            else
            {
                if (res->d == NULL)
                    res->d = new float[res->ny * res->nx];
                is >> res->d[nth_line - 2];
            }
            nth_line++;
        }
        file.close();
    }

    return res;
}

int main(int argc, char const *argv[])
{
    result *res = readFile();
    float *result = new float[res->ny * res->ny];

    correlate(res->ny, res->nx, res->d, result);
    for (int i = 0; i < res->ny * res->ny; ++i)
        cout << "i = " << i << " -- " << result[i] << endl;
    return 0;
}

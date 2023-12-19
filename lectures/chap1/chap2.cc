#include "chap2.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <random>
#include <time.h>

using namespace std;

void writeFile(string fileName, uint n)
{

    fstream file(fileName, ios::out);
    if (file.is_open())
    {
        file << n << endl;
        srand(time(0));
        for (int i = 0; i < n * n; i++)
        {
            if (i % n == i / n)
                file << 0.0 << endl;
            else
                file << static_cast<float>(rand()) / static_cast<float>(RAND_MAX) << endl;
        }

        file.close();
    }
}

chap2_result *readFile(string fileName)
{
    fstream file(fileName, ios::in);
    chap2_result *result;
    if (file.is_open())
    {
        string line;
        int i = -1;
        result = new chap2_result;
        while (getline(file, line))
        {
            istringstream iss(line);

            if (i == -1)
            {
                iss >> result->n;
                result->d = new float[result->n * result->n];
            }
            else
            {
                // cout << i << " -- " << result->d[i - 1] << endl;
                iss >> result->d[i];
            }
            ++i;
        }
        file.close();
    }

    return result;
}
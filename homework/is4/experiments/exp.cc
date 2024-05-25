#include "exp.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

// == namespace: DataGen =================================================
float DataGen::genColor()
{
    return static_cast<float>(rand() * 1.0 / RAND_MAX);
}
float *DataGen::genRandArr(int nx, int ny)
{
    float *out = new float[nx * ny * N_COLOR_CHAN];

    for (int y = 0; y < ny; ++y)
        for (int x = 0; x < nx; ++x)
            for (int c = 0; c < N_COLOR_CHAN; ++c)
            {
                out[c + 3 * x + 3 * nx * y] = DataGen::genColor();
            }

    return out;
}

void DataGen::printInput(Input &input)
{
    printf("ny = %d, nx = %d\n", input.ny, input.nx);
    for (int y = 0; y < input.ny; ++y)
        for (int x = 0; x < input.nx; ++x)
            for (int c = 0; c < N_COLOR_CHAN; ++c)
            {
                printf("%6.4f\n", input.data[c + 3 * (input.nx * y + x)]);
            }
}

Input DataGen::readFile(string filename)
{
    fstream file(filename, ios::in);
    Input input = {.ny = -1, .nx = -1, .data = nullptr};

    if (file.is_open())
    {
        string line;
        int i = -1;

        while (getline(file, line))
        {
            istringstream iss(line);

            if (i == -1)
            {
                iss >> input.ny;
                iss >> input.nx;

                input.data = new float[N_COLOR_CHAN * input.ny * input.nx];
                ++i;
            }
            else
            {
                istringstream iss(line);

                for (int j = 0; j < 3; ++j)
                {
                    iss >> input.data[i];
                    ++i;
                }
            }
        }

        file.close();
    }

    return input;
}
// == namespace: Utils =================================================

void Utils::xy2each(int xy, int nx, int &x, int &y)
{
    x = xy % nx;
    y = xy / nx;
}

int Utils::getNumCells(int left, int top, int right, int bot)
{
    int nRowsInside = bot - top + 1;
    int nColsInside = right - left + 1;

    return nRowsInside * nColsInside;
}

int Utils::getNumCellsOutside(int left, int top, int right, int bot, int nx, int ny)
{
    return Utils::getNumCells(0, 0, nx - 1, ny - 1) - Utils::getNumCells(left, top, right, bot);
}

// == namespace: Testing =================================================

template <class T>
void Testing::printColor(T *data, int nx, int ny)
{
    for (int c = 0; c < 3; ++c)
    {
        printf("- Color c = %d\n", c);

        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                printf("%8.1f", data[c + 3 * x + 3 * nx * y]);
            }
            printf("\n");
        }
    }
}
template <class T>
void Testing::printColor(T data, int nx, int ny)
{
    for (int c = 0; c < 3; ++c)
    {
        printf("- Color c = %d\n", c);

        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                printf("%8.1f", data[c + 3 * x + 3 * nx * y]);
            }
            printf("\n");
        }
    }
}

template void Testing::printColor<float>(float *data, int nx, int ny);
template void Testing::printColor<vector<double>>(vector<double> data, int nx, int ny);

void Testing::printResult(Result &result)
{
    printf("%32s\n", "OUTPUT");

    printf("%16s%16d\n", "y0", result.y0);
    printf("%16s%16d\n", "x0", result.x0);
    printf("%16s%16d\n", "y1", result.y1);
    printf("%16s%16d\n", "x1", result.x1);

    printf("%16s%+16.8f\n", "inner[0]", result.inner[0]);
    printf("%16s%+16.8f\n", "inner[1]", result.inner[1]);
    printf("%16s%+16.8f\n", "inner[2]", result.inner[2]);
    printf("%16s%+16.8f\n", "outer[0]", result.outer[0]);
    printf("%16s%+16.8f\n", "outer[1]", result.outer[1]);
    printf("%16s%+16.8f\n", "outer[2]", result.outer[2]);
}
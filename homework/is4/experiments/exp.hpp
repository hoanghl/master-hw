#ifndef EXP_HPP
#define EXP_HPP

#include <limits>
#include <math.h>
#include <string>
#include <tuple>

using namespace std;

typedef tuple<int, int, int, int> TLBR;
typedef tuple<double, double, double> COLOR;

constexpr int N_COLOR_CHAN = 3;
constexpr double inf = numeric_limits<double>::infinity();

const string FILE1 = "data/001-small-simple.txt";
const string FILE2 = "data/003-small-structured.txt";
const string FILE3 = "data/029-small-mse-map.txt";

struct Result
{
    int y0;
    int x0;
    int y1;
    int x1;
    float outer[3];
    float inner[3];
};

struct Input
{
    int ny, nx;
    float *data;
};

namespace DataGen
{
    float genColor();
    float *genRandArr(int nx, int ny);
    Input readFile(string filename);
    void printInput(Input &input);
}

namespace Utils
{
    void xy2each(int xy, int nx, int &x, int &y);
    int getNumCells(int left, int top, int right, int bot);
    int getNumCellsOutside(int left, int top, int right, int bot, int nx, int ny);
}

namespace Testing
{
    template <class T>
    void printColor(T *data, int nx, int ny);
    template <class T>
    void printColor(T data, int nx, int ny);
    void printResult(Result &result);
}

#endif
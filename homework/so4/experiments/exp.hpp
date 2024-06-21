#ifndef EXP_HPP
#define EXP_HPP

#include <string>

using namespace std;

typedef unsigned long long data_t;
typedef struct Input
{
    int n;
    data_t *data;
} Input;

const string FILE1 = "data/074-medium-threads.txt";

namespace DataGen
{
    Input *readFile(string filename);
    void genRand(int n, data_t *&data);
    void genSequential(int n, data_t *&data);
}

#endif
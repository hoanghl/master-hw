#include <iostream>
#include <string>

using namespace std;

const string STD_FILENAME = "test1.txt";

struct chap2_result
{
    int n;
    float *d;
};

#pragma once
void writeFile(string fileName, uint n);
chap2_result *readFile(string fileName);

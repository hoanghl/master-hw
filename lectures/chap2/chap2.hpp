#include <iostream>
#include <string>

using namespace std;

#ifndef CHAP2_HPP
#define CHAP2_HPP

const string STD_FILENAME = "data.txt";

struct chap2_result
{
    int n;
    float *d;
};

#pragma once
void writeFile(string fileName, uint n);
chap2_result *readFile(string fileName);

#endif
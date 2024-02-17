#include <fstream>
#include <iostream>
#include <random>
#include <string>

using namespace std;

const int N = 10000000;
const string FILENAME_FORMAT = "ex5_p1_formatted.dat";
const string FILENAME_UNFORMAT = "ex5_p1_unformatted.dat";

float getRand()
{
    return static_cast<float>(rand()) / static_cast<float>(RAND_MAX) + static_cast<float>(rand());
}

int main(int argc, char const *argv[])
{
    // Generate random numbers
    float *arr = new float[N];
    for (int i = 0; i < N; ++i)
        arr[i] = getRand();

    // Write array formatted in ASCII
    cout << "= Write array as formatted" << endl;

    clock_t start = clock();

    fstream file(FILENAME_FORMAT, ios::out);
    if (file.is_open())
    {
        for (int i = 0; i < N; ++i)
            file << arr[i];

        file.close();
    }

    clock_t end = clock();
    double elapsed = double(end - start) / CLOCKS_PER_SEC;
    cout << "CPU time (second): " << setprecision(12) << elapsed << endl;

    // Write array unformatted
    cout << "= Write array as unformatted" << endl;
    start = clock();

    fstream file2(FILENAME_UNFORMAT, ios::out | ios::binary);
    if (file2.is_open())
    {
        for (int i = 0; i < N; ++i)
            file2 << arr[i];

        file2.close();
    }

    end = clock();
    elapsed = double(end - start) / CLOCKS_PER_SEC;
    cout << "CPU time (second): " << setprecision(12) << elapsed << endl;

    return 0;
}

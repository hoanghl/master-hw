
    #include <iomanip>
    #include <iostream>
    #include <random>
    #include <time.h>

    #include "ex4.hpp"

    using namespace std;

    int main(int argc, char const *argv[])
    {
        const int N = 16777216;
        int *arr = new int[N];
        for (int i = 0; i < N; ++i)
            arr[i] = rand();
        int *sub = new int[N];

        clock_t start = clock();
        for (int i = 0; i < N; i = i + 1)
        {
    
sub[i + 0] = arr[i + 0] * 2;

}
    clock_t end = clock();
    double elapsed = double(end - start) / CLOCKS_PER_SEC;
    cout << "CPU time (second): " << setprecision(12) << elapsed << endl;

    return 0;
}

#include <iostream>

using namespace std;

int main(int argc, char const *argv[])
{
    const int N = 34;
    for (int i = 0; i < N * N; ++i)
        if (i % N == i / N)
            cout << i << endl;

    return 0;
}

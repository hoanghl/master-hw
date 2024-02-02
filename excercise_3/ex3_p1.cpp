#include <iostream>

using namespace std;

int ackerman(int m, int n)
{
    if (m == 0)
        return n + 1;
    if (m > 0 && n == 0)
        return ackerman(m - 1, 1);
    if (m > 0 && n > 0)
        return ackerman(m - 1, ackerman(m, n - 1));

    // If reach here, it has error
    cout << "Error: Either 'm' or 'n' less than zero: m = " << m << " - n = " << n << endl;
    exit(1);
}

int main()
{
    ackerman(4, 1);
}
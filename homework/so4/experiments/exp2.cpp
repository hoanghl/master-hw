#include <algorithm>
#include <iostream>
using namespace std;

// Driver code
int main()
{
    // First array
    int arr1[] = {10, 20, 30, 40, 50, 60};
    int n = sizeof(arr1) / sizeof(arr1[0]);

    // Print the first array
    cout << "Array 1: ";
    for (int i = 0; i < n; i++)
        cout << arr1[i] << " ";
    cout << endl;

    // Second array to store the elements of the first array
    int arr2[n];

    // Copy first array to second array
    copy(arr1, arr1 + n, arr2);

    // Print the second array
    cout << "Array 2: ";
    for (int i = 0; i < n; i++)
        cout << arr2[i] << " ";
    cout << endl;

    return 0;
}
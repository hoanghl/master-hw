#include <stdio.h>

int main(void)
{
    int a[10];
    int i = 1000;
    for (; i < 100000; ++i)
        printf("%d: %d\n", i, a[i]);
}
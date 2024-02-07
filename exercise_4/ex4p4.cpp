//
// Please find the example code from below.
//
// Cheers,
// Xin

#include <iostream>
#include <stdio.h>
#include <time.h>

int main()
{
  int n2 = 2560;
  int n3 = 20000000;
  int *d = new int[n2];
  int *e = new int[n2];
  int *f = new int[n2];
  clock_t v1, v2;
  /*---------Initialize vectors---------*/
  v1 = clock();
  for (int i = 0; i < n2; i++)
  {
    d[i] = 5;
    e[i] = 5 + 3;
    f[i] = 0;
  }
  v2 = clock();
  printf("Initialize_vetor: \t %.6f s\n", (double)(v2 - v1) / CLOCKS_PER_SEC);
  clock_t v3, v4;
  /*---------f = d * e (repeat)---------*/
  v3 = clock();
  for (int i = 0; i < n3; i++)
  {
    for (int j = 0; j < n2; j++)
    {
      f[j] = d[j] * e[j];
    }
  }
  v4 = clock();
  printf("de(repeat): \t\t\t %.6f s\n", (double)(v4 - v3) / CLOCKS_PER_SEC);
  delete[] d;
  delete[] e;
  delete[] f;
  return 0;
}

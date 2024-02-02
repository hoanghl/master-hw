#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define N 15000
int a[N][N];
#define NPRINT 500

void main()
{
  int i, j;
  clock_t t1, t2, tcpu;

  t1 = clock();
  /* Begin measurement */
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      a[i][j] = (i + j) / 2;
  /* End measurement */
  t1 = clock();
  tcpu = (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
  printf("CPU time: %f\n", (double)tcpu);

  for (i = 0; i < N; i += NPRINT)
    for (j = 0; j < N; j += NPRINT)
      printf("%d ", a[i][j]);
  printf("\n");
}

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NMAX 100000
#define MMAX 1000

void main() {
  double dx=0.00001;
  double sum1,sum2;
  int i,j;

  sum1=sum2=0.0;
  
  for (j=0;j<=MMAX;j++) 
    for (i=0;i<=NMAX;i++) 
      sum1=sum1+dx*i*j;

  for (j=0;j<=MMAX;j++) 
    for (i=0;i<=NMAX;i++) 
      sum2=sum2+cos(exp(sin(dx*i*j)));

  printf("%g %g\n",sum1,sum2);

}

#include <stdlib.h>
#include <stdio.h>

#define N 15000
int a[N][N];
#define NPRINT 500

void main() {
  int i,j;
  clock_t t1,t2;


  /* Begin measurement */
  for (i=0;i<N;i++) 
    for (j=0;j<N;j++) 
      a[i][j]=(i+j)/2;
  /* End measurement */


  for (i=0;i<N;i+=NPRINT)
    for (j=0;j<N;j+=NPRINT)
      printf("%d ",a[i][j]);
  printf("\n");


}


  
  
    

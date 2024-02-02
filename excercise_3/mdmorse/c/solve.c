/* $Id: solve.c 18 2008-09-16 09:54:33Z aakurone $ */

#include "global.h"

void Solve1(struct vector  *x,
	    struct vector  *v,
	    struct vector  *a,
	    int             N,
	    struct vector   size,
	    struct ivector  periodic,
	    double          deltat)

{
  double   halfdeltatsq;
  struct vector   half;

  int  i;
  
  halfdeltatsq=0.5*deltat*deltat;

  half.x=size.x/2.0;
  half.y=size.y/2.0;
  half.z=size.z/2.0;

  for(i=1;i<=N;i++) {

     /* First part of velocity Verlet solution */
     
     x[i].x=x[i].x+deltat*v[i].x+halfdeltatsq*a[i].x;
     x[i].y=x[i].y+deltat*v[i].y+halfdeltatsq*a[i].y;
     x[i].z=x[i].z+deltat*v[i].z+halfdeltatsq*a[i].z;

     v[i].x=v[i].x+0.5*deltat*a[i].x;
     v[i].y=v[i].y+0.5*deltat*a[i].y;
     v[i].z=v[i].z+0.5*deltat*a[i].z;
     
     /* Enforce periodics if needed */

     if (periodic.x==1 && x[i].x < -half.x) x[i].x=x[i].x+size.x;
     if (periodic.x==1 && x[i].x >= half.x) x[i].x=x[i].x-size.x;
     
     if (periodic.y==1 && x[i].y < -half.y) x[i].y=x[i].y+size.y;
     if (periodic.y==1 && x[i].y >= half.y) x[i].y=x[i].y-size.y;
     
     if (periodic.z==1 && x[i].z < -half.z) x[i].z=x[i].z+size.z;
     if (periodic.z==1 && x[i].z >= half.z) x[i].z=x[i].z-size.z;

  }

}


void Solve2(struct vector  *v,
	    struct vector  *a,
	    int             N,
	    double          deltat,
	    double          btctau,
	    double          T,
	    double          desiredT)
  
{
  int  i;
  double btclambda;

  /* 
     Compute btc lambda. See J. Chem. Phys. 81 (1984) 3684.
     Note the factor 2 inside the sqrt. For details see the lecture notes.
  */

  if (btctau > 0.0 && T > 0.0) {
     btclambda = sqrt(1+2.0*deltat/btctau*(desiredT/T-1.0));
  } else {
     btclambda=1.0;
  }
  
  for(i=1;i<=N;i++) {

     /* Second part of velocity Verlet solution */
     
     v[i].x=v[i].x+0.5*deltat*a[i].x;
     v[i].y=v[i].y+0.5*deltat*a[i].y;
     v[i].z=v[i].z+0.5*deltat*a[i].z;

     v[i].x=v[i].x*btclambda;
     v[i].y=v[i].y*btclambda;
     v[i].z=v[i].z*btclambda;

  }

}


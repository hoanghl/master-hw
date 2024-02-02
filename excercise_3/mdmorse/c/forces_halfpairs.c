/* $Id: forces_halfpairs.c 18 2008-09-16 09:54:33Z aakurone $ */

#include "global.h"
  
void GetForces(struct vector   *x,
	       struct vector   *a,
	       int              N,
	       struct vector    size,
	       struct ivector   periodic,
	       double           m,
	       int             *neighbourlist,
	       double           rpotcut, 
	       double           morseD,
	       double           morsealpha,
	       double           morser0,
	       double          *Epot,
	       double          *virial)

{
  struct vector   half;
  int        i,j,jj,startofineighbourlist,nneighboursi;
  double   r,rsq,rpotcutsq,dx,dy,dz,help1,help2,V,dVdr;
  
  half.x=size.x/2.0;
  half.y=size.y/2.0;
  half.z=size.z/2.0;
  rpotcutsq=rpotcut*rpotcut;

  for(i=1;i<=N;i++) {
     Epot[i]=0.0;
     a[i].x=0.0;
     a[i].y=0.0;
     a[i].z=0.0;
  }
  
  *virial=0.0;
  startofineighbourlist=1;
  for(i=1;i<=N-1;i++) {
     nneighboursi=neighbourlist[startofineighbourlist];

     /* Loop over neighbours of i in neighbourlist */

     for (jj=1;jj<=nneighboursi;jj++) {
       /* Retrive real atom index */
       j=neighbourlist[startofineighbourlist+jj];

       /* Get distance */
       dx=x[j].x-x[i].x;
       if (periodic.x==1 && dx>half.x ) dx=dx-size.x;
       if (periodic.x==1 && dx<-half.x) dx=dx+size.x;

       dy=x[j].y-x[i].y;
       if (periodic.y==1 && dy>half.y ) dy=dy-size.y;
       if (periodic.y==1 && dy<-half.y) dy=dy+size.y;

       dz=x[j].z-x[i].z;
       if (periodic.z==1 && dz>half.z ) dz=dz-size.z;
       if (periodic.z==1 && dz<-half.z) dz=dz+size.z;

       rsq=dx*dx+dy*dy+dz*dz;

       if (rsq <= rpotcutsq) {

	   /* i-j interaction within cutoff, evaluate it */
  
	   r=sqrt(rsq);

	   /* Utilize the fact that one of the Morse terms is the */
	   /* square of the other */
	   help1=exp(-morsealpha*(r-morser0));
	   help2=help1*help1;

	   /* Calculate potential energy */
	   V=0.5*morseD*(help2-2.0*help1);
	   Epot[i]=Epot[i]+V;
	   Epot[j]=Epot[j]+V;

	   /* Calculate forces. Take care with signs */

	   /* First get dV/dr in units of eV/Å */
	   dVdr=-morseD*(help2*(-2.0*morsealpha)-2.0*help1*(-morsealpha));

	   /* Get virial = Sum (r·f) for pressure calculation in units of eV  */
	   *virial=*virial+dVdr*(dx/r*dx+dy/r*dy+dz/r*dz);

	   /* then do transformations to get dV/dr in SI units */
	   dVdr=dVdr*e/lunit;
	   /* Then get a=F/m in units of Å/fs^2 */
	   dVdr=(dVdr/(m*u))/aunit/r;

	   /* Now project on each vector */

	   a[i].x=a[i].x-dVdr*dx;
	   a[j].x=a[j].x+dVdr*dx;
	   a[i].y=a[i].y-dVdr*dy;
	   a[j].y=a[j].y+dVdr*dy;
	   a[i].z=a[i].z-dVdr*dz;
	   a[j].z=a[j].z+dVdr*dz;
	   
       }

     } /* End of loop over neighbours jj */

     /* Set starting position for next atom */
     startofineighbourlist=startofineighbourlist+nneighboursi+1;

  } /* End of loop over atoms i */


}

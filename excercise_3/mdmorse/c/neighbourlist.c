/* $Id: neighbourlist.c 18 2008-09-16 09:54:33Z aakurone $ */

#include "global.h"

void UpdateNeighbourlist(struct vector   *x,
			 int              N,
			 struct vector    size,
			 struct ivector   periodic,
			 int             *neighbourlist,
			 double           rskincut)


/* 

   Subroutine constructs a Verlet neighbour list in the following format:

   nngbr1          Number of neighbours for atom 1
   ingbr1,1        Index of neighbour 1 of atom 1
   ingbr2,1        Index of neighbour 2 of atom 1
   ...
   ingbrnngbr1,1   Index of neighbour nngbr1 of atom 1
   nngbr2          Number of neighbours for atom 2
   
*/

/* NOTE: The first index is 1, not 0 !! */
{

  double   dx,dy,dz,rsq,rskincutsq;
  struct vector   half;
  int      i,j,nneighboursi,startofineighbourlist,nneighbourstot;
  
  half.x=size.x/2.0;
  half.y=size.y/2.0;
  half.z=size.z/2.0;
  rskincutsq=rskincut*rskincut;;

  nneighbourstot=0;
  startofineighbourlist=1;
  for(i=1;i<=N;i++) {
     nneighboursi=0;
     
     for(j=1;j<=N;j++) {	
	if (i==j) continue;
	
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

	if (rsq <= rskincutsq) {
	   nneighboursi=nneighboursi+1;
	   nneighbourstot=nneighbourstot+1;
	   if (startofineighbourlist+nneighboursi > MAXAT*MAXNEIGHBOURS) {
	     puts("Neighbour list overflow, increase MAXNEIGHBOURS");
	     exit(0);
	   }
	   neighbourlist[startofineighbourlist+nneighboursi]=j;
	}
		
     } /* End of loop over neighbours j */

     /* Write in number of i's neighbours into list */
     neighbourlist[startofineighbourlist]=nneighboursi;
     /* Set starting position for next atom */
     startofineighbourlist=startofineighbourlist+nneighboursi+1;
     
  } /* End of loop over atoms i */

  printf("%s %13.5f %s\n","Neighbour list update found",
       1.0*nneighbourstot/N,"neighbours per atom");
  
}



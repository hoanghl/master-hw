/* $Id: main.c 20 2008-09-16 13:17:57Z aakurone $ */


/* A simple MD program made for demonstration purposes.

   Internal units:
   Length:       Å  = 1e-10 m
   Times:        fs = 1e-15 s
   Velocity:     Å/fs
   Accelration:  Å/fs^2
   Forces:       eV/Å
   Energy:       eV

   NOTE: Since this code is converted from Fortran, array
   indices run from 1 to N, not from 0 to N-1 as usually in C. 

*/

#include "global.h"

main()
{  

  /* 

     Arrays containing atom data:
     
     x,v,a    : atom positions, velocities and accelerations
     Epot     : potential energy for each atom
     Ekin     : kinetic energy     - " - 
     atomname : simply the atom name 

*/

  struct vector   x[MAXAT+1],v[MAXAT+1],a[MAXAT+1];
  double   Epot[MAXAT+1],Ekin[MAXAT+1];
  char atomname[MAXAT+1][5];
  int      neighbourlist[MAXAT*MAXNEIGHBOURS+1];
  
  /* 

     Other parameters and variables
  
     time                 : current time
     tmax                 : maximum time
     deltat               : time step
     itime                : integer number of the current time step
     T                    : temperature (in K)
     m                    : atom mass
     Ekinsum, Epotsum     : system kinetic and potential energies
     Etot                 : total energy for all atoms
     N                    : number of atoms
     nneighbourlistupdate : interval between neighbour list updates
     nmovieoutput         : interval between atom position output
     rpotcut              : potential cutoff radius
     rskincut             : neighbour list cutoff radius
     size                 : simulation cell size
     periodic             : periodic boundaries will be used (1=yes, 0=no) 
     seed                 : random number generator seed

*/
  
  double   time,tmax,deltat,T,Ekinsum,Epotsum,Etot,m;
  int      N,nneighbourlistupdate,itime,nmovieoutput;
  double   rpotcut,rskincut;
  struct vector   size;
  struct ivector  periodic;
  int seed;

  /* Morse potential parameters */
  double  morseD,morsealpha,morser0;

  /*

    Parameters for Berendsen temperature and pressure control
    [J. Chem. Phys.81 (1984) 3684]

    btctau   : Berendsen temperature control tau,
    initialT : initial temperature
    desiredT : desired temperature
    bpctau   : Berendsen pressure control tau
    bpcbeta  : 1/B
    P        : instantaneous pressure
    desiredP : desired pressure

  */

  double btctau,initialT,desiredT;
  double bpctau,bpcbeta,P,desiredP,virial,mu;


  int  i;
  
  printf("\n--------------- mdmorse08 1.0c --------------------\n\n");

  /* Initialize the simulation */

  ReadParams(&initialT,&m,&tmax,&deltat,&size,&periodic,
	     &nneighbourlistupdate,&rpotcut,&rskincut,  
	     &morseD,&morsealpha,&morser0,&nmovieoutput,
	     &btctau,&desiredT,&bpctau,&bpcbeta,&desiredP,&seed);
   
  ReadAtoms(x,atomname,&N);
  
  /* Factor 2 because of energy equipartition theorem */
  SetTemperature(v,N,m,2.0*initialT,&seed);

  /* Check that we got it right */
  GetTemperature(v,N,m,&T);
  printf("Initial atom temperature is %g\n",T);

  /* Initialize accelerations to 0 */
  for (i=1;i<=N;i++) {
     a[i].x=0.0; a[i].y=0.0; a[i].z=0.0;
  }
  
  /* Start main loop over times */
  time=0.0;
  itime=0;
  while (1) {
    if (itime%nneighbourlistupdate==0) {
      UpdateNeighbourlist(x,N,size,periodic,neighbourlist,rskincut);
    }
     
    Solve1(x,v,a,N,size,periodic,deltat);

    GetForces(x,a,N,size,periodic,m,neighbourlist,rpotcut, 
	      morseD,morsealpha,morser0,Epot,&virial);

    GetTemperature(v,N,m,&T);

    Solve2(v,a,N,deltat,btctau,T,desiredT);

    /* Update time */
    time=time+deltat;
     
    /* Actual solution done, now compute result quantities */

    GetEnergies(v,N,m,&Ekinsum,&Epotsum,&Etot,Ekin,Epot);

    /* Compute pressure in units of kbar */
    P=(N*kB*T+1.0/3.0*virial*e)/(size.x*size.y*size.z*1e-30)/1e8;

    printf("%s %10.3f %g %12.5f %12.5f %12.5f %12.5f\n","ec ",
	   time,T,Ekinsum/N,Epotsum/N,Etot/N,P);

    if (nmovieoutput>0) 
      if (itime%nmovieoutput==0) {
	printf("%s %10.3f\n","Outputting atom movie at t = ",time);
	WriteAtoms(x,atomname,Ekin,Epot,N,time,size,itime);
      }
    
    /* Do Berendsen pressure control */
    if (bpctau > 0.0) {
      mu=pow((1.0-bpcbeta*deltat/bpctau*(desiredP-P)),(1.0/3.0));

      size.x=size.x*mu;
      size.y=size.y*mu;
      size.z=size.z*mu;
      for (i=1;i<=N;i++) {
	x[i].x=x[i].x*mu;
	x[i].y=x[i].y*mu;
	x[i].z=x[i].z*mu;
      }
      if (itime%10==1) {
	printf("bpc %10.3f %10.6f %10.6f %10.6f %12.5f %12.5f %12.5f\n",
	       time,size.x,size.y,size.z,size.x*size.y*size.z,P,mu);
      }
    }

    /* Check whether we are done */
    if (time >= tmax) break;

    itime=itime+1;

  } /* End of loop over time steps */

  /* Always output final energy and atom coordinates using high precision. */
  printf("%s %12.5f %12.5f %14.7f %14.7f %14.7f\n","ec ",
	 time,T,Ekinsum/N,Epotsum/N,Etot/N);

  WriteAtoms(x,atomname,Ekin,Epot,N,time,size,itime);

}
  

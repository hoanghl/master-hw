/* $Id: physical.c 20 2008-09-16 13:17:57Z aakurone $ */


#include "global.h"

void SetTemperature(struct vector  *v,
		    int             N,
		    double          m,
		    double          T,
		    int            *seed)

{  


  /* 

     Subroutine sets the temperature of the atoms in the system
     to T assuming a Maxwell-Boltzmann velocity distribution. 

  */


  struct vector  vtot;
  double  vxrms,gaussianrand();
  int  i;
  
  /* RMS v in one dimension. See e.g. Mandl p. 191 */
  
  vxrms=sqrt(kB*T/(m*u))/vunit;

  vtot.x=0.0; vtot.y=0.0; vtot.z=0.0; 
  for (i=1;i<=N;i++) {

     v[i].x=vxrms*gaussianrand(seed);
     v[i].y=vxrms*gaussianrand(seed);
     v[i].z=vxrms*gaussianrand(seed);

     vtot.x=vtot.x+v[i].x;
     vtot.y=vtot.y+v[i].y;
     vtot.z=vtot.z+v[i].z;
     
  }

  vtot.x=vtot.x/N;
  vtot.y=vtot.y/N;
  vtot.z=vtot.z/N;

  /* Subtract total velocity of atoms  */
  for (i=1;i<=N;i++) {
     v[i].x=v[i].x-vtot.x;
     v[i].y=v[i].y-vtot.y;
     v[i].z=v[i].z-vtot.z;
  }
  
}


void GetTemperature(struct vector  *v,
		    int             N,
		    double          m,
		    double         *T)


{
  /*

    Subroutine gets the temperature of the atoms in the system.

  */

  double  Ekinsum,gaussianrand();
  int  i;

  /* Get sum of kinetic energies */
  Ekinsum=0;
  for (i=1;i<=N;i++) {

     Ekinsum=Ekinsum+v[i].x*v[i].x+v[i].y*v[i].y+v[i].z*v[i].z;

  }
  /* Do necessary unit transformations to get E in J */
  Ekinsum=0.5*m*u*Ekinsum*vunit*vunit;
  /* and get temperature using E=3/2 N k T */
  *T = Ekinsum/(1.5*N*kB);
  
}

void GetEnergies(struct vector  *v,
		 int             N,
		 double          m,
		 double         *Ekinsum,
		 double         *Epotsum,
		 double         *Etot,
		 double         *Ekin,
		 double         *Epot)

{
  double  help1;
  int  i;
  
  *Ekinsum=0.0; *Epotsum=0.0;
  help1=0.5*m*u*vunit*vunit/e;
  for (i=1;i<=N;i++) {
     Ekin[i]=help1*(v[i].x*v[i].x+v[i].y*v[i].y+v[i].z*v[i].z);
     *Ekinsum = *Ekinsum+Ekin[i];
     *Epotsum = *Epotsum+Epot[i];
  }
  *Etot = *Ekinsum+*Epotsum;

}
  



/*--------------------------------------------------------------------------- */

double uniformrand(int *seed)

  /*

    Park-Miller "minimal" Random number generator uniform distribution ]0,1[. 
    See Numerical Recipes ch. 7.0.

  */

{
  int IA,IM,IQ,IR,MASK;
  double AM,ans;
  int  k;
  
  IA=16807; IM=2147483647; AM=1.0/IM; IQ=127773; IR=2836; MASK=123459876;

  *seed= (*seed)^MASK;
  k=(*seed)/IQ;
  *seed = IA*((*seed)-k*IQ)-IR*k;
  if (*seed < 0) *seed = *seed+IM;
  ans=AM*(*seed);
  *seed= (*seed)^MASK;
  
  return(ans);
}



double gaussianrand(int *seed)

  /*

    Random numbers with normal (Gaussian) distribution.
    Mean is 0 and standard deviation is 1.
    See W.H.Press et al., Numerical Recipes 1st ed., page 203

  */

{
  double fac,v1,v2,r,ans;
  static double gset;
  static int iset=0;
  double uniformrand();

  if (iset==0) {
    do {
      v1 = 2.0*uniformrand(seed)-1.0;
      v2 = 2.0*uniformrand(seed)-1.0;
      r = v1*v1+v2*v2;
    } while (r>=1.0||r==0.0);
    fac = sqrt(-2.0*log(r)/r);
    gset = v1*fac;
    iset = 1;
    return(v2*fac); 
  }
  else {
    iset = 0;
    return(gset);
  }

}

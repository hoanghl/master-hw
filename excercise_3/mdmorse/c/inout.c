/* $Id: inout.c 22 2008-11-12 13:06:32Z aakurone $ */

/* Input-output subroutines for the mdmorse code */

#include "global.h"

#define MAXBUF 200

void ReadParams(double         *initialT,
		double         *m,
		double         *tmax,
		double         *deltat,
		struct vector  *size,
		struct ivector *periodic,
		int            *nneighbourlistupdate,
		double         *rpotcut,
		double         *rskincut,  
		double         *morseD,
		double         *morsealpha,
		double         *morser0,
		int            *nmovieoutput,
		double         *btctau,
		double         *desiredT,
		double         *bpctau,
		double         *bpcbeta,
		double         *desiredP,
		int            *seed)


{
  int  i,n;
  FILE *fp;
  char  string[MAXBUF],buf[MAXBUF],*c;
  double x;
  
  fp=fopen("mdmorse.in","r");
   
  if (fp==NULL) { puts("Open error for file mdmorse.in"); exit(1); }

  /* Initialize variables to sensible (?) default values */
  *initialT=0.0; *m=63.546;
  size->x=0.0; size->y=0.0; size->z=0.0;

  *deltat=1.0;
  *nneighbourlistupdate=5;
  *rpotcut=4.02; *rskincut=5.0;
  /* Initialize Morse parameters for Cu */
  /* Values from Girifalco and Weizer, Phys. Rev. 114 (1959) 687. */
  *morseD=0.3429; *morsealpha=1.3588; *morser0=2.866;
  *btctau=0.0; *desiredT=0.0;
  *bpctau=0.0; *bpcbeta=1e-3; *desiredP=0.0;

  *seed=23267761;


  
  /* Loop through file and find all variables */
  i=0;
  while(1) {
     i=i+1;
     c=fgets(buf,MAXBUF,fp);
     if (c==NULL) break;

     string[0]=' ';
     n=sscanf(buf,"%s %lg",string,&x);

     /* If not parameter line, skip this line */
     if (string[0] != '$') continue;

     if (n<2) {
       printf("ERROR: Line \"%s\" looks like data line\n",buf);
       printf("but does not contain data in correct format\n");
       exit(0);
     }

     /* Look for variable-defining strings */
     if (strcmp(string,"$initialT")==0) *initialT=x;
     else if (strcmp(string,"$desiredT")==0) *desiredT=x;
     else if (strcmp(string,"$btctau")==0) *btctau=x;

     else if (strcmp(string,"$mass")==0) *m=x;

     else if (strcmp(string,"$bpctau")==0) *bpctau=x;
     else if (strcmp(string,"$bpcbeta")==0) *bpcbeta=x;
     else if (strcmp(string,"$desiredP")==0) *desiredP=x;

     else if (strcmp(string,"$xsize")==0) size->x=x;
     else if (strcmp(string,"$ysize")==0) size->y=x;
     else if (strcmp(string,"$zsize")==0) size->z=x;

     else if (strcmp(string,"$periodicx")==0) periodic->x=(int)(x+0.5);
     else if (strcmp(string,"$periodicy")==0) periodic->y=(int)(x+0.5);
     else if (strcmp(string,"$periodicz")==0) periodic->z=(int)(x+0.5);
     
     else if (strcmp(string,"$seed")==0) *seed=(int)(x+0.5);
     
     else if (strcmp(string,"$tmax")==0) *tmax=x;
     else if (strcmp(string,"$deltat")==0) *deltat=x;

     else if (strcmp(string,"$nupdate")==0) *nneighbourlistupdate=(int)(x+0.5);

     else if (strcmp(string,"$rpotcut")==0) *rpotcut=x;
     else if (strcmp(string,"$rskincut")==0) *rskincut=x;
     
     else if (strcmp(string,"$morseD")==0) *morseD=x;
     else if (strcmp(string,"$morsealpha")==0) *morsealpha=x;
     else if (strcmp(string,"$morser0")==0) *morser0=x;

     else if (strcmp(string,"$nmovieoutput")==0) *nmovieoutput=(int)(x+0.5);

     else {
       printf("Unknown parameter %s in input file\n",string);
       exit(0);
     }

     printf("Read in parameter %s to value %g\n",string,x);

  } /* End of readin loop */

  fclose(fp);

  /* Print out some general information about read results: */

  printf("%s %3d %3d %3d\n","Using periodics (1=on, 0=off)",
	 periodic->x,periodic->y,periodic->z);

  printf("%s %12.6g %12.6g %12.6g\n","Morse potential parameters: D alpha r0",
       *morseD,*morsealpha,*morser0);

  printf("%s %d %s\n","Movie output selected every",*nmovieoutput," steps");

  if (*btctau > 0.0) {
    printf("Doing Berendsen temperature control with tau T %g %g\n",
	   *btctau,*desiredT);
  }

  if (*bpctau > 0.0) {
    printf("Doing Berendsen pressure control with tau beta %g %g\n",
	   *bpctau,*bpcbeta);
  }

  return;
  
}



void ReadAtoms(struct vector *x,
	       char           atomname[MAXAT+1][5],
	       int           *N)

{
  char comment[MAXBUF],buf[MAXBUF];
  int  i,n;
  FILE *fp;

  fp=fopen("atoms.in","r");

  if (fp==NULL) { puts("File atoms.in open failed"); exit(0); }

  fgets(buf,MAXBUF,fp); sscanf(buf,"%d",N);

  if (*N > MAXAT) { puts("Error: MAXAT overflow, increase it."); exit(0); }

  fgets(comment,MAXBUF,fp);
  printf("%s%d%s\"%s\"\n","Reading in ",*N," atoms described as ",comment);
  
  for (i=1;i<= *N;i++) {

    fgets(buf,MAXBUF,fp); 
    n=sscanf(buf,"%s %lg %lg %lg",atomname[i],
	   &(x[i].x),&(x[i].y),&(x[i].z));

    if (n<4) {
	printf("Atom read in error for line \"%s\" n %d\n",buf,n);
	printf("%s %lg %lg %lg\n",atomname[i],x[i].x,x[i].y,x[i].z);
	puts("Atom read in error"); 
	exit(0);
    }

  }
  
  
}



void WriteAtoms(struct vector  *x,
		char            atomname[MAXAT+1][5],
		double         *Ekin,
		double         *Epot,
		int             N,
		double          time,
		struct vector   size,
		int             itime)
  
{
  
  char  comment[MAXBUF];
  int  i;
  static FILE *fp;

  static int firsttime=1;

  if (firsttime) {
    fp=fopen("atoms.out","w");
    if (fp==NULL) { puts("File atoms.out open failed"); exit(0); }
    firsttime=0;
  }

  fprintf(fp," %d\n",N);
  
  /*
  fprintf(fp,"%s %11.3f %s %12.4f %12.4f %12.4f\n",
	  "mdmorse atom output at time",
  	  time,"fs boxsize ",size.x,size.y,size.z);
  */
  fprintf(fp," Frame number %d %g fs boxsize %g %g %g\n",itime,time,size.x,size.y,size.z);
  
  for(i=1;i<=N;i++) {
    fprintf(fp,"%2s %12.6f%12.6f%12.6f%11.6f\n",
	    atomname[i],x[i].x,x[i].y,x[i].z,Epot[i]);
  }

  fflush(fp);
  
}

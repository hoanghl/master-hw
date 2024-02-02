/* $Id: global.h 20 2008-09-16 13:17:57Z aakurone $ */


/* Define a vector type */
struct vector {
  double  x,y,z;
};

/* Define an int vector type */
struct ivector {
  int  x,y,z;
};

#define  pi 3.14159265358979
#define  e 1.6021892e-19
#define  u 1.660565e-27
#define  kB 1.380662e-23

/* 
   
   Define some unit conversion factors.
   Needed e.g. in setting and calculating the temperature.
   To get internal units from SI, divide by one of these.
   To get SI from internal units, multiply by one of these.

*/

#define  lunit 1e-10
#define  tunit 1e-15
#define  vunit (lunit/tunit)
#define  aunit (vunit/tunit)
  
/* Set array sizes */

#define  MAXAT 10000
#define  MAXNEIGHBOURS 100
  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

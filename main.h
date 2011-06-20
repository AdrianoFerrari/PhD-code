#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <mpi.h>
#include "bessels.h"

#define PI 3.14159265358979323846
#define twoPI 6.283185307179586

/* KISS RNG 
https://groups.google.com/group/comp.lang.fortran/browse_thread/thread/a85bf5f2a97f5a55?fwc=2&hl=en
*/

#define MWC (t=(x<<58)+c, c=(x>>6), x+=t, c+=(x<t), x) 
#define XSH ( y^=(y<<13), y^=(y>>17), y^=(y<<43) ) 
#define CNG ( z=6906969069LL*z+1234567 ) 
#define KISS (MWC+XSH+CNG) 
#define UNI (KISS*5.421010862e-20)
typedef unsigned long long ULL;
static ULL x=1234567890987654321ULL,c=123456123456123456ULL, y=362436362436362436ULL,z=1066149217761810ULL,t; 

void settable(ULL i1, ULL i2, ULL i3, ULL i4)
{ 
  x = i1; c = i2; y = i3; z = i4;
}

const double uy = 0.041666666666;
const double Ly = 24.0;
const int M = 7;//Lekner terms

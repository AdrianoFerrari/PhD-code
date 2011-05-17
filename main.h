#include <math.h>
#include "bessels.h"

#define PI 3.14159265358979323846
#define twoPI 6.283185307179586

/* KISS RNG 
https://groups.google.com/group/comp.lang.fortran/browse_thread/thread/a85bf5f2a97f5a55?fwc=2&hl=en
*/
static unsigned long long 
x=1234567890987654321ULL,c=123456123456123456ULL, 
  y=362436362436362436ULL,z=1066149217761810ULL,t; 
#define MWC (t=(x<<58)+c, c=(x>>6), x+=t, c+=(x<t), x) 
#define XSH ( y^=(y<<13), y^=(y>>17), y^=(y<<43) ) 
#define CNG ( z=6906969069LL*z+1234567 ) 
#define KISS (MWC+XSH+CNG) 

const double uy = 0.041666666666;
const int M = 7;//Lekner terms

//double dummy_function();
//double lekner_potential(double q0, double x0, double y0, double z0, double q1, double x1, double y1, double z1 );
//double repulsive_potential(double x0,double y0,double z0,double x1,double y1,double z1, double epsilon);

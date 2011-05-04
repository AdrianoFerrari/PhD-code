#include <math.h>
#include "bessels.h"

#define PI 3.14159265358979323846
#define twoPI 6.283185307179586

const double uy = 0.041666666666;
const int M = 7;//Lekner terms

double dummy_function();
double lekner_potential(double q0, double x0, double y0, double z0, double q1, double x1, double y1, double z1 );

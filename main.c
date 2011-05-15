#include "main.h"

double dummy_function() {
  return 1.234;
}
double lekner_potential(double q0, double x0, double y0, double z0, double q1, double x1, double y1, double z1 ) {
  double rxz = sqrt((x1 - x0) * (x1 - x0) + (z1 - z0) * (z1 - z0));
  double y = (y1 - y0);
  double vi = 0;
  int n;
  for (n = 1; n <= M; n++)
    {
      vi += 4.0*q0*q1*cos(twoPI*n*y*uy)*bessk0(twoPI*n*rxz*uy)*uy;
    }
  return vi - 2.0*q0*q1*log(rxz)*uy;
}
double repulsive_potential(double x0,double y0,double z0,double x1,double y1,double z1, double epsilon) {
  double r2 = (x1 - x0)*(x1 - x0) + (y1-y0)*(y1-y0) + (z1 - z0)*(z1 - z0);
  if (r2 > 10.0*pow(epsilon,0.166666)) {
    return 0.0;
  }
  else if (r2 == 0.0) {
    return 10e20;
  }
  else {
    return epsilon/pow(r2,6);
  }
}

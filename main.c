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
double line_potential(double q1, double x1, double y1, double z1, double qL, double xL, double zL, double A, double lambda) {
  double rxz = sqrt((x1 - xL) * (x1 - xL) + (z1 - zL) * (z1 - zL));
  if(rxz<=0)
    return 0.0;
  else
    return -2.0*q1*qL*log(rxz)+
		4.0*A*PI*q1*x1*qL*bessk1(twoPI*rxz/lambda)*sin(twoPI*y1/lambda)/rxz;
}
double line_repulsive_potential(double epsilon, double x0,double y0,double z0,double xL,double zL, double A, double lambda) {
  double xa = xL + A*sin(twoPI*y0/lambda);
  double r2 = (x0 - xa)*(x0 - xa) + (z0 - zL)*(z0 - zL);
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
double total_line_potential(double x[][3], double q[], double eps, int i, double qL, double xL, double zL, double A, double lambda) {
  return line_potential(q[i],x[i][0],x[i][1],x[i][2],qL,xL,zL,A,lambda) + line_repulsive_potential(eps,x[i][0],x[i][1],x[i][2],xL,zL,A,lambda);
}
double pair_potential_energy(double q0, double x0, double y0, double z0, double q1, double x1, double y1, double z1, double epsilon) {
  return repulsive_potential(x0,y0,z0,x1,y1,z1,epsilon) + lekner_potential(q0,x0,y0,z0,q1,x1,y1,z1);
}
double particle_total_pair_potential(double x[][3], double q[], int size, double eps, int i) {
  int j;
  double pe = 0.0;
  double q0 = q[i];
  double x0 = x[i][0];  
  double y0 = x[i][1];  
  double z0 = x[i][2];

  for(j=0;j<size;j++) {
    if(j != i) {
      pe += pair_potential_energy(q0,x0,y0,z0,q[j],x[j][0],x[j][1],x[j][2],eps);
    }
  }
  
  return pe;
}
double particle_total_potential(double x[][3], double q[], int size, int i, double eps, double qL, double xL, double zL, double A, double lambda) {
  return particle_total_pair_potential(x,q,size,eps,i)+total_line_potential(x,q,eps,i,qL,xL,zL,A,lambda);
}
double energy_difference(double x[][3], double xn[][3], double q[], int size, int i, double eps, double qL, double xL, double zL, double A, double lambda) {
  double energy = particle_total_potential(x,q,size,i,eps,qL,xL,zL,A,lambda);
  double new_energy = particle_total_potential(xn,q,size,i,eps,qL,xL,zL,A,lambda);
  return new_energy - energy;
}

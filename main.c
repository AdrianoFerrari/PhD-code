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
double particle_total_potential(double x[][3], double q[], int size, int i, double eps, double qL, double xL1, double zL1, double A1, double xL2, double zL2, double A2, double lambda) {
  return particle_total_pair_potential(x,q,size,eps,i)+total_line_potential(x,q,eps,i,qL,xL1,zL1,A1,lambda)+total_line_potential(x,q,eps,i,qL,xL2,zL2,A2,lambda);
}
double energy_difference(double x[][3], double xn[][3], double q[], int size, int i, double eps, double qL, double xL1, double zL1, double A1, double xL2, double zL2, double A2, double lambda) {
  double energy = particle_total_potential(x,q,size,i,eps,qL,xL1,zL1,A1,xL2,zL2,A2,lambda);
  double new_energy = particle_total_potential(xn,q,size,i,eps,qL,xL1,zL1,A1,xL2,zL2,A2,lambda);
  return new_energy - energy;
}
bool go_ahead(double De, double kbT) {
  if(UNI < fmin(1.0,exp(-De/kbT))) {
      return true;
    }
  else return false;
}
double line_growth_at_pt_from_particle(double s, double q, double x, double y, double z, double qL, double xL, double zL, double A, double lambda, int Ns) {
  double dx = (x - xL - A*sin(twoPI*s*Ly/lambda));
  double rxz = sqrt(dx*dx  + (z - zL) * (z - zL));
  double G = 0.0;
  double vi = 0.0;
  int n = 1;
  for (n = 1; n <= 12; n++)
    {
      if (rxz > 0) vi += 8.0*PI*qL*q*n*cos(twoPI*n*(y-s*Ly)*uy)*bessk1(twoPI*n*rxz*uy)*uy/(1.0*Ns);
    }
  G = vi + 2.0*qL*q/(rxz*Ns);
  return G/(A*sin(twoPI*s*Ly/lambda)*sqrt(1+(z-zL)*(z-zL)/(dx*dx)));
}
double xforce_on_seg_due_to_part(double s,double q,double x,double y,double z,double qL,double xL,double zL,double ep, double A,double lambda, int Ns,int M) {
  double xs = xL + A*sin(twoPI*s*Ly/lambda);
  double rxz = sqrt((x-xs)*(x-xs)+(z-zL)*(z-zL));
  double vi = 0.0;
  int n = 1;
  for (n =1; n <= M; n++) {
      if (rxz > 0) vi += 8.0*PI*qL*q*n*cos(twoPI*n*(y-s*Ly)*uy)*bessk1(twoPI*n*rxz*uy)*uy/(1.0*Ns);
    }
  double frxz = vi + 12*ep*pow(rxz,-13) + 2.0*qL*q/(rxz*Ns);
  if (x-xs >= 0.0) {
    return -1.0*frxz/sqrt(1+(z-zL)*(z-zL)/((x-xs)*(x-xs)));
  }
  else {
    return frxz/sqrt(1+(z-zL)*(z-zL)/((x-xs)*(x-xs)));
  } 
}
double xforce_on_seg_due_to_line(double s, double qL, double xL1, double zL1, double xL2, double zL2, double A, double lambda, int Ns) {
  double xs = xL1 + A*sin(twoPI*s*Ly/lambda);
  double rxz = sqrt((xL2-xs)*(xL2-xs)+(zL2-zL1)*(zL2-zL1));
  double frxz = 2*qL*qL*Ly/(rxz*Ns);
  if (xL2-xs >= 0.0) {
    return -1.0*frxz/sqrt(1+(zL2-zL1)*(zL2-zL1)/((xL2-xs)*(xL2-xs)));
  }
  else {
    return frxz/sqrt(1+(zL2-zL1)*(zL2-zL1)/((xL2-xs)*(xL2-xs)));
  }
}
double growth_on_seg_due_to_part(double s, double q, double x, double y, double z, double qL, double xL, double zL, double ep, double A, double lambda, int Ns, int M) {
  return (xforce_on_seg_due_to_part(s,q,x,y,z,qL,xL,zL,ep,A,lambda,Ns,M)-xforce_on_seg_due_to_part(s,q,x,y,z,qL,xL,zL,ep,0.0,lambda,Ns,M))/(A*sin(twoPI*s*Ly/lambda));
}
double growth_on_seg_due_to_line(double s, double qL, double xL1, double zL1, double xL2, double zL2, double A, double lambda, int Ns) {
  return (xforce_on_seg_due_to_line(s,qL,xL1,zL1,xL2,zL2,A,lambda,Ns)-xforce_on_seg_due_to_line(s,qL,xL1,zL1,xL2,zL2,0.0,lambda,Ns))/(A*sin(twoPI*s*Ly/lambda));
}

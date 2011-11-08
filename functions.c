#include "main.h"

double dist(y1,y0) {
  if abs(y1-y0) > Ly/2.0 {
    return Ly/2.0 - abs(y1-y0);
  }
  else {
    return abs(y1-y0);
  }
}

double lekner_u(double q0, double x0, double y0, double z0, double q1, double x1, double y1, double z1 ) {
  double rxz = sqrt((x1 - x0) * (x1 - x0) + (z1 - z0) * (z1 - z0));
  double y = dist(y1,y0);
  double vi = 0;
  int n;
  for (n = 1; n <= M; n++)
    {
      if (rxz != 0) {
	vi += 4.0*q0*q1*cos(twoPI*n*y*uy)*bessk0(twoPI*n*rxz*uy)*uy;}
    }
  return vi - 2.0*q0*q1*log(rxz)*uy;
}
double repulsive_u(double ep, double x0,double y0,double z0,double x1,double y1,double z1) {
  double r2 = (x1 - x0)*(x1 - x0) + dist(y1,y0)*dist(y1,y0) + (z1 - z0)*(z1 - z0);
  if (r2 > 10.0*pow(epsilon,0.166666)) {
    return 0.0;
  }
  else if (r2 == 0.0) {
    return 10e30;
  }
  else {
    return epsilon/pow(r2,6);
  }
}
double spring_u(double h,double Lmax,double x0,double y0,double z0,double x1,double y1,double z1) {
  L = sqrt((x1 - x0)*(x1 - x0) + dist(y1,y0)*dist(y1,y0) + (z1 - z0)*(z1 - z0));
  if (L >= Lmax) {
    return 1e10;
  }
  else {
    return -(0.5*h*Lmax*Lmax)*log(1-(L/Lmax)*(L/Lmax));
  }
}
bool on_chain(int i,int N,int Nl){
  if (i < N) { return false; }
  else {return true; }
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
double xforce_on_seg(double s, double q[], double x[][3], int size,double qL, double xL1, double zL1, double xL2, double zL2, double ep, double A, double lambda, int Ns, int M) {
  double fLine = xforce_on_seg_due_to_line(s,qL,xL1,zL1,xL2,zL2,A,lambda,Ns);
  double fParts = 0.0;
  int i = 0;
  for(i = 0; i < size; i++) {
    fParts += xforce_on_seg_due_to_part(s,q[i],x[i][0],x[i][1],x[i][2],qL,xL1,zL1,ep,A,lambda,Ns,M);
  }
  return fLine + fParts;
}
double xforce(double q[], double x[][3], int size,double qL, double xL1, double zL1, double xL2, double zL2, double ep, double A, double lambda, int Ns, int M) {
  double sum = 0.0;
  int i = 1;
  for(i = 1;i <= Ns; i++) {
    sum += xforce_on_seg(i*Ly/(1.0*Ns),q,x,size,qL,xL1,zL1,xL2,zL2,ep,A,lambda,Ns,M);
  }
  return sum;
}
double growth_on_seg(double s, double q[], double x[][3], int size,double qL, double xL1, double zL1, double xL2, double zL2, double ep, double A, double lambda, int Ns, int M) {
  return (xforce_on_seg(s,q,x,size,qL,xL1,zL1,xL2,zL2,ep,A,lambda,Ns,M)-xforce_on_seg(s,q,x,size,qL,xL1,zL1,xL2,zL2,ep,0.0,lambda,Ns,M))/(A*sin(twoPI*s*Ly/lambda));
}
double growth_rate(double q[], double x[][3], int size,double qL, double xL1, double zL1, double xL2, double zL2, double ep, double A, double lambda, int Ns, int M) {
  double sum = 0.0;
  int i = 1;
  for(i = 1; i <=Ns; i++) {
    sum += growth_on_seg(i*Ly/(1.0*Ns),q,x,size,qL,xL1,zL1,xL2,zL2,ep,A,lambda,Ns,M);
  }
  return sum/(1.0*Ns);
}
double ran_xz(double maxR) { return sqrt(maxR*maxR*UNI)*cos(twoPI*UNI); }
double ran_y(double maxY) { return maxY*UNI; }
int ran_particle(int total) { return (int)floor(total*UNI); }

#include "main.h"

double lekner_u(double q0, double x0, double y0, double z0, double q1, double x1, double y1, double z1 ) {
  double rxz = sqrt((x1 - x0) * (x1 - x0) + (z1 - z0) * (z1 - z0));
  double y = dist(y1,y0);
  double vi = 0;
  int n;
  for (n = 1; n <= M; n++)
    {
      if (rxz != 0) {
	vi += 4.0*q0*q1*cos(twoPI*n*y*uy)*gsl_sf_bessel_K0(twoPI*n*rxz*uy)*uy;}
    }
  return vi - 2.0*q0*q1*log(rxz)*uy;
}
double repulsive_u(double ep, double x0,double y0,double z0,double x1,double y1,double z1) {
  double r2 = (x1 - x0)*(x1 - x0) + dist(y1,y0)*dist(y1,y0) + (z1 - z0)*(z1 - z0);
  if (r2 > 10.0*pow(ep,0.166666)) {
    return 0.0;
  }
  else if (r2 == 0.0) {
    return 10e30;
  }
  else {
    return ep/pow(r2,6);
  }
}
double spring_u(double h,double Lmax,double x0,double y0,double z0,double x1,double y1,double z1) {
  double L = sqrt((x1 - x0)*(x1 - x0) + dist(y1,y0)*dist(y1,y0) + (z1 - z0)*(z1 - z0));
  if (L >= Lmax) {
    return 1e10;
  }
  else {
    return -(0.5*h*Lmax*Lmax)*log(1-(L/Lmax)*(L/Lmax));
  }
}
bool on_chain(int i,int N,int Nl) {
  if (i < N) { return false; }
  else {return true; }
}
bool is_endpoint(int i,int N,int Nl) {
  if (i == N || i == N+Nl-1 || i == N+2*Nl-1) {
    return true;
  }
  else {return false;}
}
/*bool linked(int i, int j, int N, int Nl) {
  if (i>j) {
    int t;
    t = i;
    i = j;
    j = t;
  }
  if (!on_chain(i,N,Nl) || !on_chain(j,N,Nl)) {
    return false;
  }
  else if ((i==N && j==N+Nl-1)||(i==N+Nl && j == N+2*Nl-1)) {
    return true;
  }
  else if (abs(i-j) != 1) { return false; }
  else { return true; }
  }*/
double total_u(int i,double q0,double x0,double y0,double z0,double j,double q1,double x1,double y1,double z1,double ep,double h,double Lmax,int N,int Nl) {
  double total = 0.0;
  total += lekner_u(q0,x0,y0,z0,q1,x1,y1,z1) + repulsive_u(ep,x0,y0,z0,x1,y1,z1);
  
  if(linked(i,j,N,Nl)) {
    total += spring_u(h,Lmax,x0,y0,z0,x1,y1,z1);
  }
  
  return total;
}
double delta_u(double **x, double **xn, double *q, int n, double ep, double h, double Lmax, int N, int Nl) {
  double delta_e = 0.0;
  
#pragma omp parallel shared(n,N,Nl,Lmax,ep,h,q,x,xn) reduction(+:delta_e)
  {
    #pragma omp for
    for (int i = 0; i < N+2*Nl; i++)
      {
	if(i!=n) {
	  double new_e, old_e;
	  double lek = 0; double rep = 0; double spring = 0;
	  double rxz, r, y;
	
      //-OLD config energy
	  rxz = sqrt((x[i][0]-x[n][0])*(x[i][0]-x[n][0]) + (x[i][2]-x[n][2])*(x[i][2]-x[n][2]));
	  y   = dist(x[i][1],x[n][1]);
	  r   = sqrt(rxz*rxz+y*y);
	  
	  //---Lekner E
	  for(int n = 1; n <=M; n++)
	    {
	      if(rxz != 0) {
		lek += 4.0*q[n]*q[i]*cos(twoPI*n*y*uy)*gsl_sf_bessel_K0(twoPI*n*rxz*uy)*uy;}
	    }
	  lek = lek - 2.0*q[n]*q[i]*log(rxz)*uy;
	  
	  //---Repulsive E
	  if(r*r <= 10.0*pow(ep,0.1666666) && r != 0) {
	    rep = ep/pow(r*r,6);
	  }
	  else if (r == 0.0) {
	    rep = 10e30;
	  }

	  //---Spring E
	  if(linked(i,n,N,Nl) || linked(n,i,N,Nl)) {
	    if (r >= Lmax) { spring = 1e10; }
	    else { spring = -(0.5*h*Lmax*Lmax)*log(1-(r/Lmax)*(r/Lmax)); }
	  }
	  
	  //---Sum old E
	  old_e = lek + rep + spring;

      //-NEW config energy
	  lek = rep = spring = 0;
	  rxz = sqrt((xn[i][0]-xn[n][0])*(xn[i][0]-xn[n][0]) + (xn[i][2]-xn[n][2])*(xn[i][2]-xn[n][2]));
	  y   = dist(xn[i][1],xn[n][1]);
	  r   = sqrt(rxz*rxz+y*y);
	  
	  //---Lekner E
	  for(int n = 1; n <=M; n++)
	    {
	      if(rxz != 0) {
		lek += 4.0*q[n]*q[i]*cos(twoPI*n*y*uy)*gsl_sf_bessel_K0(twoPI*n*rxz*uy)*uy;}
	    }
	  lek = lek - 2.0*q[n]*q[i]*log(rxz)*uy;
	  
	  //---Repulsive E
	  if(r*r <= 10.0*pow(ep,0.1666666) && r != 0) {
	    rep = ep/pow(r*r,6);
	  }
	  else if (r == 0.0) {
	    rep = 10e30;
	  }

	  //---Spring E
	  if(linked(i,n,N,Nl) || linked(n,i,N,Nl)) {
	    if (r >= Lmax) { spring = 1e10; }
	    else { spring = -(0.5*h*Lmax*Lmax)*log(1-(r/Lmax)*(r/Lmax)); }
	  }
	  
	  //---Sum old E
	  new_e = lek + rep + spring;
	  
     //-DELTA E
	  delta_e += new_e - old_e;
	}
      }
  }
  return delta_e;
}
double ran_xz(double maxR) { return sqrt(maxR*maxR*UNI)*cos(twoPI*UNI); }
double ran_y(double maxY) { return maxY*UNI; }
int ran_particle(int total) { return (int)floor(total*UNI); }
double ran_u() { return UNI; }
double ran_du() { return 2.0*UNI - 1.0; }

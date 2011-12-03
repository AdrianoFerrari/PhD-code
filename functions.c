#include "main.h"

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
	    if (r >= Lmax) { spring = 10e30; }
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
	    if (r >= Lmax) { spring = 10e30; }
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

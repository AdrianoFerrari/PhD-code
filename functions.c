#include "main.h"

bool on_chain(int i,int N,int Nl) {
  if (i < N) { return false; }
  else {return true; }
}
bool is_endpoint(int i,int N,int Nl) {
  if (i == N || i == N+Nl-1 || i == N+Nl || i == N+2*Nl-1) {
    return true;
  }
  else {return false;}
}
double theta_du(double x1, double y1, double z1, double x0, double y0, double z0, double xp, double yp, double zp, double xn, double yn, double zn){
  double ao = sqrt( (x0-xp)*(x0-xp) + dist(y0,yp)*dist(y0,yp) + (z0-zp)*(z0-zp) );
  double bo = sqrt( (x0-xn)*(x0-xn) + dist(y0,yn)*dist(y0,yn) + (z0-zn)*(z0-zn) );
  double c2 = (xn-xp)*(xn-xp) + dist(yn,yp)*dist(yn,yp) + (zn-zp)*(zn-zp);
  double an = sqrt( (x1-xp)*(x1-xp) + dist(y1,yp)*dist(y1,yp) + (z1-zp)*(z1-zp) );
  double bn = sqrt( (x1-xn)*(x1-xn) + dist(y1,yn)*dist(y1,yn) + (z1-zn)*(z1-zn) );

  double cos_th0 = (ao*ao+bo*bo-c2)/(2.0*ao*bo);
  double cos_th1 = (an*an+bn*bn-c2)/(2.0*an*bn);
  double th0     =  cos_th0 <= -1.0 ? PI
                  : cos_th0 >=  1.0 ? 0.0
                  : acos(cos_th0);
  double th1     =  cos_th1 <= -1.0 ? PI
                  : cos_th1 >=  1.0 ? 0.0
                  : acos(cos_th1);
  
  return th1*th1-th0*th0-twoPI*(th1-th0);
}
double delta_u(double **x, double **xn, double *q, int n, double ep, double h, double htheta, double Lmax, int N, int Nl) {
  double delta_e = 0.0;
  double x1 = xn[n][0];
  double y1 = xn[n][1];
  double z1 = xn[n][2];
  double x0 = x[n][0];
  double y0 = x[n][1];
  double z0 = x[n][2];

#pragma omp parallel shared(n,N,Nl,Lmax,ep,h,q,x,xn,x1,y1,z1,x0,y0,z0) reduction(+:delta_e)
  {
    #pragma omp for
    for (int i = 0; i < N+2*Nl; i++)
      {
	if(i!=n) {
	  double new_e, old_e;
	  double lek = 0; double rep = 0; double spring = 0;
	  double rxz, r, y;
	
      //-OLD config energy
	  rxz = sqrt((x[i][0]-x0)*(x[i][0]-x0) + (x[i][2]-z0)*(x[i][2]-z0));
	  y   = dist(x[i][1],y0);
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
	  rxz = sqrt((xn[i][0]-x1)*(xn[i][0]-x1) + (xn[i][2]-z1)*(xn[i][2]-z1));
	  y   = dist(xn[i][1],y1);
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
  
  //---Delta Spring Theta
  double dtheta;

  if(n>=N) {
    if(n==N){
      dtheta = theta_du(x1,y1,z1,x0,y0,z0,x[N+Nl-1][0],x[N+Nl-1][1],x[N+Nl-1][2],x[n+1][0],x[n+1][1],x[n+1][2]);
    } else if(n==N+Nl-1) {
      dtheta = theta_du(x1,y1,z1,x0,y0,z0,x[n-1][0],x[n-1][1],x[n-1][2],x[N][0],x[N][1],x[N][2]);
    } else if(n==N+Nl) {
      dtheta = theta_du(x1,y1,z1,x0,y0,z0,x[N+2*Nl-1][0],x[N+2*Nl-1][1],x[N+2*Nl-1][2],x[n+1][0],x[n+1][1],x[n+1][2]);
    } else if(n==N+2*Nl-1) {
      dtheta = theta_du(x1,y1,z1,x0,y0,z0,x[n-1][0],x[n-1][1],x[n-1][2],x[N+Nl][0],x[N+Nl][1],x[N+Nl][2]);
    } else {
      dtheta = theta_du(x1,y1,z1,x0,y0,z0,x[n-1][0],x[n-1][1],x[n-1][2],x[n+1][0],x[n+1][1],x[n+1][2]);
    }
    delta_e += htheta*dtheta;
  }
  
  return delta_e;
}
double ran_xz(double maxR) { return sqrt(maxR*maxR*UNI)*cos(twoPI*UNI); }
double ran_y(double maxY) { return maxY*UNI; }
int ran_particle(int total) { return (int)floor(total*UNI); }
double ran_u() { return UNI; }
double ran_du() { return 2.0*UNI - 1.0; }

// Force calculation functions

// force on particle 0 by particle 1
double lekner_fx(double q0, double x0, double y0, double z0, double q1, double x1, double y1, double z1) {
  double fxz, x, rxz, y, r;
  fxz = 0.0;
  x   = (x0-x1);
  rxz = sqrt(x*x + (z0-z1)*(z0-z1));
  y   = dist(y0,y1);
  r   = sqrt(rxz*rxz+y*y);

  for(int n=1; n<=M; n++) {
    if(rxz == 0) continue;
    fxz += 4*twoPI*q1*q0*n*cos(twoPI*n*y*uy)*gsl_sf_bessel_K1(twoPI*n*rxz*uy)*uy*uy;
  }
  fxz = fxz + 2.0*q1*q0*uy/rxz;

  return fxz*x/r;
}

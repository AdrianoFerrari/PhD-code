#include "main.h"

double dist(double y1, double y0, double Ly) {
  return fabs((y1-y0)-Ly*round((y1-y0)/Ly));
}
bool on_chain(int i,int N) {
  if (i < N) { return false; }
  else {return true; }
}
bool is_endpoint(int i,int N,int Nl) {
  if (i == N || i == N+Nl-1 || i == N+Nl || i == N+2*Nl-1) {
    return true;
  }
  else {return false;}
}
double energy(double **x, double *q, double ep, double sigma, double sigma_c, double h, double htheta, double Lmax, int N, int Nl, double Ly, Bilin_interp lekner) {
  double totalE = 0.0;

  for (int i = 0; i < N+2*Nl; i++) {
    for (int j = i+1; j < N+2*Nl; j++) {
      double lek = 0; double rep = 0; double spring = 0;
      double uy = 1.0/Ly;
      double rxz, r, y;
    
      rxz = sqrt((x[i][0]-x[j][0])*(x[i][0]-x[j][0]) + (x[i][2]-x[j][2])*(x[i][2]-x[j][2]));
      y   = dist(x[i][1],x[j][1],Ly);
      r   = sqrt(rxz*rxz+y*y);
      
      lek = 2.0*q[i]*q[j]*uy*(lekner.interp(rxz*uy,y*uy)-log(Ly));
      
      if(j>=N && i>=N) {
        if(r <= 1.12246*sigma_c) {
          double sigma_r_6 = pow(sigma_c/r,6);
          rep = 4.0*ep*(sigma_r_6*sigma_r_6-sigma_r_6+0.25);
        }
        else {
          rep = 0.0;
        }
      }
      else {
        if(r <= 1.12246*sigma) {
          double sigma_r_6 = pow(sigma/r,6);
          rep = 4.0*ep*(sigma_r_6*sigma_r_6-sigma_r_6+0.25);
        }
        else {
          rep = 0.0;
        }
      }

      if(linked(i,j,N,Nl) || linked(j,i,N,Nl)) {
        if (r >= Lmax) { spring = INFINITY; }
        else { spring = -0.5*h*Lmax*Lmax*log(1-r*r/(Lmax*Lmax)); }
      }
    
    totalE += lek + rep + spring;
    }
  }

  return totalE;
}
double theta_du(double x1, double y1, double z1, double x0, double y0, double z0, double xp, double yp, double zp, double xn, double yn, double zn, double Ly){
  double ao = sqrt( (x0-xp)*(x0-xp) + dist(y0,yp,Ly)*dist(y0,yp,Ly) + (z0-zp)*(z0-zp) );
  double bo = sqrt( (x0-xn)*(x0-xn) + dist(y0,yn,Ly)*dist(y0,yn,Ly) + (z0-zn)*(z0-zn) );
  double c2 = (xn-xp)*(xn-xp) + dist(yn,yp,Ly)*dist(yn,yp,Ly) + (zn-zp)*(zn-zp);
  double an = sqrt( (x1-xp)*(x1-xp) + dist(y1,yp,Ly)*dist(y1,yp,Ly) + (z1-zp)*(z1-zp) );
  double bn = sqrt( (x1-xn)*(x1-xn) + dist(y1,yn,Ly)*dist(y1,yn,Ly) + (z1-zn)*(z1-zn) );

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
double delta_u(double **x, double **xn, double *q, int n, double ep, double sigma, double sigma_c, double h, double htheta, double Lmax, int N, int Nl, double Ly, Bilin_interp lekner) {
  double x1 = xn[n][0];
  double y1 = xn[n][1];
  double z1 = xn[n][2];
  double x0 = x[n][0];
  double y0 = x[n][1];
  double z0 = x[n][2];

  double uy = 1.0/Ly;
  double new_e, old_e, lek, rep, spring, rxz, y, r;

  double delta_e = 0.0;

#pragma omp parallel for private(new_e, old_e, lek, rep, spring, rxz, y, r) reduction(+:delta_e)
  for (int i = 0; i < N+2*Nl; i++) {
    if(i!=n) {
      double new_e, old_e;
      double lek = 0; double rep = 0; double spring = 0;
      double rxz, y, r;
    
        //-OLD config energy
      rxz = sqrt((x[i][0]-x0)*(x[i][0]-x0) + (x[i][2]-z0)*(x[i][2]-z0));
      y   = dist(x[i][1],y0,Ly);
      r   = sqrt(rxz*rxz+y*y);
      
      //---Lekner E
      lek = 2.0*q[n]*q[i]*uy*(lekner.interp(rxz*uy,y*uy)-log(Ly));
      
      //---Repulsive E
      if(n>=N && i>=N) {
        if(r <= 1.12246*sigma_c) {
          double sigma_r_6 = pow(sigma_c/r,6);
          rep = 4.0*ep*(sigma_r_6*sigma_r_6-sigma_r_6+0.25);
        }
        else {
          rep = 0.0;
        }
      }
      else {
        if(r <= 1.12246*sigma) {
          double sigma_r_6 = pow(sigma/r,6);
          rep = 4.0*ep*(sigma_r_6*sigma_r_6-sigma_r_6+0.25);
        }
        else {
          rep = 0.0;
        }
      }

      //---Spring E
      if(linked(i,n,N,Nl) || linked(n,i,N,Nl)) {
        if (r >= Lmax) { spring = INFINITY; }
        else { spring = -0.5*h*Lmax*Lmax*log(1-r*r/(Lmax*Lmax)); }
      }
      
      //---Sum old E
      old_e = lek + rep + spring;

        //-NEW config energy
      lek = rep = spring = 0;
      rxz = sqrt((xn[i][0]-x1)*(xn[i][0]-x1) + (xn[i][2]-z1)*(xn[i][2]-z1));
      y   = dist(xn[i][1],y1,Ly);
      r   = sqrt(rxz*rxz+y*y);
      
      //---Lekner E
      lek = 2.0*q[n]*q[i]*uy*(lekner.interp(rxz*uy,y*uy)-log(Ly));
      
      //---Repulsive E
      if(n>=N && i>=N) {
        if(r <= 1.12246*sigma_c) {
          double sigma_r_6 = pow(sigma_c/r,6);
          rep = 4.0*ep*(sigma_r_6*sigma_r_6-sigma_r_6+0.25);
        }
        else {
          rep = 0.0;
        }
      }
      else {
        if(r <= 1.12246*sigma) {
          double sigma_r_6 = pow(sigma/r,6);
          rep = 4.0*ep*(sigma_r_6*sigma_r_6-sigma_r_6+0.25);
        }
        else {
          rep = 0.0;
        }
      }

      //---Spring E
      if(linked(i,n,N,Nl) || linked(n,i,N,Nl)) {
        if (r >= Lmax) { spring = INFINITY; }
        else { spring = -(0.5*h*Lmax*Lmax)*log(1-(r/Lmax)*(r/Lmax)); }
      }
      
      //---Sum old E
      new_e = lek + rep + spring;
      
      //-DELTA E
      delta_e += new_e - old_e;
    }
  }

  //---Delta Spring Theta
  double dtheta;

  if(n>=N) {
    if(n==N){
      dtheta = theta_du(x1,y1,z1,x0,y0,z0,x[N+Nl-1][0],x[N+Nl-1][1],x[N+Nl-1][2],x[n+1][0],x[n+1][1],x[n+1][2],Ly);
    } else if(n==N+Nl-1) {
      dtheta = theta_du(x1,y1,z1,x0,y0,z0,x[n-1][0],x[n-1][1],x[n-1][2],x[N][0],x[N][1],x[N][2],Ly);
    } else if(n==N+Nl) {
      dtheta = theta_du(x1,y1,z1,x0,y0,z0,x[N+2*Nl-1][0],x[N+2*Nl-1][1],x[N+2*Nl-1][2],x[n+1][0],x[n+1][1],x[n+1][2],Ly);
    } else if(n==N+2*Nl-1) {
      dtheta = theta_du(x1,y1,z1,x0,y0,z0,x[n-1][0],x[n-1][1],x[n-1][2],x[N+Nl][0],x[N+Nl][1],x[N+Nl][2],Ly);
    } else {
      dtheta = theta_du(x1,y1,z1,x0,y0,z0,x[n-1][0],x[n-1][1],x[n-1][2],x[n+1][0],x[n+1][1],x[n+1][2],Ly);
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

// Lekner force_x on particle 0 by particle 1
double lekner_fx(double q0, double x0, double y0, double z0, double q1, double x1, double y1, double z1, double Ly) {
  double fxz, x, rxz, y, r;
  double uy = 1.0/Ly;
  fxz = 0.0;
  x   = (x0-x1);
  rxz = sqrt(x*x + (z0-z1)*(z0-z1));
  y   = dist(y0,y1,Ly);
  r   = sqrt(rxz*rxz+y*y);
  if(rxz == 0.0 && y != 0.0) { return 0.0;}
  if(rxz == 0.0 && y == 0.0) { return NAN;}

  for(int n=1; n<=M; n++) {
    fxz += 4*twoPI*q1*q0*n*cos(twoPI*n*y*uy)*gsl_sf_bessel_K1(twoPI*n*rxz*uy)*uy*uy;
  }
  fxz = fxz + 2.0*q1*q0*uy/rxz;

  return fxz*x/r;
}
// Replusive force_x on particle 0 by particle 1
double rep_fx(double ep, double x0, double y0, double z0, double x1, double y1, double z1, double Ly) {
  double x = (x0-x1);
  double r = sqrt( x*x + dist(y0,y1,Ly)*dist(y0,y1,Ly) + (z0-z1)*(z0-z1) );
  if(r*r <= 10.0*pow(ep,0.1666666) && r != 0) {
    return 12.0*ep*x/pow(r,14);
  }
  else {
    return 0.0;
  }
}
// Spring force_x on particle 0 by particle 1
double spring_fx(double h, double Lmax, double x0, double y0, double z0, double x1, double y1, double z1, double Ly) {
  double x = (x0-x1);
  double r = sqrt( x*x + dist(y0,y1,Ly)*dist(y0,y1,Ly) + (z0-z1)*(z0-z1) );
  if(r < Lmax) { return -1.0*h*x/(1-r*r/(Lmax*Lmax)); }
  else { return 0.0; }
}

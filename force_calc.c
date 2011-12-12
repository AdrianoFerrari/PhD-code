0;136;0c#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_bessel.h>

#define PI 3.14159265358979323846
#define twoPI 6.283185307179586

const double Ly = 24.0;
const double uy = 0.041666666666;
const int M = 7;//Lekner terms

#define dist(y1,y0) Ly*fabs(fabs(y1-y0)/Ly-floor(fabs(y1-y0)/Ly+0.5))

double lekner_fx(double q0, double x0, double y0, double z0, double q1, double x1, double y1, double z1) {
  double fxz, x, rxz, y, r;
  fxz = 0.0;
  x   = (x1-x0);
  rxz = sqrt(x*x + (z1-z0)*(z1-z0));
  y   = dist(y1,y0);
  r   = sqrt(rxz*rxz+y*y);

  for(int n=1; n<=M; n++) {
    if(rxz == 0) continue;
    fxz += 4*twoPI*q1*q0*n*cos(twoPI*n*y*uy)*gsl_sf_bessel_K1(twoPI*n*rxz*uy)*uy*uy;
  }
  fxz = fxz + 2.0*q1*q0*uy/rxz;

  return fxz*x/r;
}
												       

int main(int argc, char **argv) {
  if(argc != 6) { printf("Usage: N Nl ci ep h\n"); }

  int N             = atoi(argv[1]);
  int Nl            = atoi(argv[2]);
  double ci_charge  = atof(argv[3]);
  double ep         = atof(argv[4]);
  double h          = atof(argv[5]);

  //inits
  double qL = -0.5*N*ci_charge;

  //init arrays
  double *q;
  double **x;
  q = malloc((N+2*Nl)*sizeof(double));
  x = malloc((N+2*Nl)*sizeof(double *));
  
  for(int i=0; i< N+2*Nl; i++) {
    x[i] = malloc(3*sizeof(double));
    if(x[i] == NULL) { printf("Out of memory!\n"); }
  }

  //open input file
  FILE *input;
  input = fopen("coords.xyz","r");
  char buff[1024];
  int ni; double xc,yc,zc;

  int i = 0;
  while(fgets(buff, 1024, input)) {
    sscanf(buff,"%d %lf %lf %lf", &ni, &xc, &yc, &zc);
    x[i][0] = xc;
    x[i][1] = yc;
    x[i][2] = zc;
    i++;
  }
  

  //initialize charges
  for(int i=0; i<N+2*Nl;i++) {
    if(i<N) {
      q[i] = ci_charge;
    }
    else {
      q[i] = qL/(1.0*Nl);
    }
  }

  //caclulate lekner
  double fLx = 0.0;
  double fRx = 0.0;
  for(int i=0; i<N+2*Nl; i++) {
    for(int j=0; j< N+2*Nl;j++) {
      if(i==j) { continue; }
      
      if(i >= N && i < N+Nl) {//if i is on left chain, add to fLx
	fLx += lekner_fx(q[i],x[i][0],x[i][1],x[i][2],q[j],x[j][0],x[j][1],x[j][2]);
      }
      else if(i >= N+Nl) {//if i is on the right chain, add to fRx
	fRx += lekner_fx(q[i],x[i][0],x[i][1],x[i][2],q[j],x[j][0],x[j][1],x[j][2]);
      }
    }
  }

  //open output file
  FILE *output;
  output = fopen("force.xyzf","w");

  fprintf(output,"fLx\tfRx\n");
  fprintf(output,"%f\t%f\n",fLx,fRx);
	
  //cleanup & close
  fclose(output);
  for(int i=0; i< N+2*Nl; i++) {
    free(x[i]);
  }
  free(x); free(q);
}

#include "functions.c"
// #include <stdio.h>

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

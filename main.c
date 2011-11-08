#include "tests.c"

int main(int argc, char **argv) {
  run_tests(); 
  
  //variable init
  int s;
  double kbt;
  char filename[32];
  int i;  double De; int ranN;

  //import simulation parameters
  int N             = atoi(argv[1]);
  int Nl            = atoi(argv[2]);
  double ci_charge  = atof(argv[3]);
  double ep         = atof(argv[4]);
  int T             = atoi(argv[4]);
  double kbt0       = atof(argv[5]);
  double kbtf       = atof(argv[6]);
  double R          = atof(argv[7]);
  sprintf(filename,"%s",argv[8]);
  int posOut        = atoi(argv[9]);
  
  //constants & inits
  double maxR = 16.0;
  double maxY = 24.0;
  double qL = -0.5*N*ci_charge;

  //init arrays
  double q[N];
  double x[N][3];
  double xn[N][3];

  //set random seed
  settable(1234567890987654321ULL, 123456123456123456ULL, 362436362436362436ULL, 1066149217761810ULL);

  //open required output files
  char fname [32];
  FILE *pos;
  if (posOut != 0) {
    sprintf(fname, "%s.xyz",filename); pos=fopen(fname,"w");
  }

  //initialize particles
  for(i=0; i<N;i++) {
    q[i] = ci_charge;
    x[i][0] = ran_xz(maxR);
    x[i][1] = ran_y(maxY);
    x[i][2] = ran_xz(maxR);
    xn[i][0] = x[i][0];
    xn[i][1] = x[i][1];
    xn[i][2] = x[i][2];
  }

  //MC loop
  for(s=0; s < T; s++) {
    ranN = ran_particle(N);
    xn[ranN][0] = ran_xz(maxR);
    xn[ranN][1] = ran_y(maxY);
    xn[ranN][2] = ran_xz(maxR);
    
    kbt = kbt0;
    
    if(go_ahead(energy_difference(x,xn,q,N,ranN,ep,qL,), kbtf)) {
      x[ranN][0] = xn[ranN][0];
      x[ranN][1] = xn[ranN][1];
      x[ranN][2] = xn[ranN][2];
    }
    else {
      xn[ranN][0] = x[ranN][0];
      xn[ranN][1] = x[ranN][1];
      xn[ranN][2] = x[ranN][2];
    }
    
    if(posOut != 0 && s % posOut == 0) {
      fprintf(pos,"%d\n",N);
      fprintf(pos,"somestuff\n");
      for(i=0;i<N;i++) {
	fprintf(pos,"0 %f %f %f\n",x[i][0],x[i][1],x[i][2]);
      }
    }
  }

  //close files
  if(posOut != 0) { fclose(pos); }
}

#include "functions.c"

int main(int argc, char **argv) {
  if(argc != 13) { printf("Usage: ./main N Nl qci ep h Lmax T kbt0 kbtf R file posOut\n"); return 0; }
  
  //variable init
  int s;
  double kbt;
  double dx, dy, dz, step;
  double acceptance_rate = 0.0;
  double De;
  char filename[32];
  int i; int ranN;

  //import simulation parameters
  int N             = atoi(argv[1]);
  int Nl            = atoi(argv[2]);
  double ci_charge  = atof(argv[3]);
  double ep         = atof(argv[4]);
  double h          = atof(argv[5]);
  double Lmax       = atof(argv[6]);
  int T             = atoi(argv[7]);
  double kbt0       = atof(argv[8]);
  double kbtf       = atof(argv[9]);
  double R          = atof(argv[10]);
  sprintf(filename,"%s",argv[11]);
  int posOut        = atoi(argv[12]);
  
  //inits
  double qL = -0.5*N*ci_charge;

  //init arrays
  double q[N+2*Nl];
  double x[N+2*Nl][3];
  double xn[N+2*Nl][3];

  //set random seed
  settable(1234567890987654321ULL, 123456123456123456ULL, 362436362436362436ULL, 1066149217761810ULL);

  //open required output files
  char fname [32];
  FILE *pos;
  if (posOut != 0) {
    sprintf(fname, "%s.xyz",filename); pos=fopen(fname,"w");
  }

  //initialize particles
  for(i=0; i<N+2*Nl;i++) {
    if(i<N) {
      q[i] = ci_charge;
      x[i][0] = xn[i][0] = ran_xz(maxR);
      x[i][1] = xn[i][1] = ran_y(Ly);
      x[i][2] = xn[i][2] = ran_xz(maxR);
    }
    else if (i<N+Nl) {
      q[i] = qL;
      x[i][0] = xn[i][0] = -0.5*R;
      x[i][1] = xn[i][1] = (i-N)*Ly/Nl - 0.5*Ly;
      x[i][2] = xn[i][2] = 0.0;
    }
    else {
      q[i] = qL;
      x[i][0] = xn[i][0] = 0.5*R;
      x[i][1] = xn[i][1] =  (i-N-Nl)*Ly/Nl - 0.5*Ly;
      x[i][2] = xn[i][2] = 0.0;
    }
  }

  //MC loop
  for(s=0; s < T; s++) {
    kbt = kbt0;
    ranN = ran_particle(N+2*Nl);
    dx = ran_u(); dy = ran_u(); dz = ran_u();
    step = on_chain(ranN, N, Nl) ? 0.2 : 0.6;
    xn[ranN][0] += step*dx;
    xn[ranN][1] += is_endpoint(ranN,N,Nl) ? 0.0 : step*dy;
    xn[ranN][2] += step*dz;
    
    De = delta_u(xn,x,q,ranN,ep,h,Lmax,N,Nl);
    
    if(De < 0 || exp(-De/kbt) > ran_u()) {
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
      fprintf(pos,"%d\n",N+2*Nl);
      fprintf(pos,"rundata\n");
      for(i=0;i<N+2*Nl;i++) {
	fprintf(pos,"%d %f %f %f\n",on_chain(i,N,Nl)?10:1,x[i][0],x[i][1],x[i][2]);
      }
    }
  }

  //close files
  if(posOut != 0) { fclose(pos); }
}

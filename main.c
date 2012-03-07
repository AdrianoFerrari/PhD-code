#include "functions.c"

int main(int argc, char **argv) {
  if(argc != 15) { printf("Usage: ./main N Nl qci ep h htheta Lmax T kbt R file posOut forceOut seed\n"); return 0; }
  pca_time tt;
  tick(&tt);
  
  //variable init
  int s;
  double dx, dy, dz, step, fLx, fRx;
  double fLavg = 0.0;
  double fRavg = 0.0;
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
  double htheta     = atof(argv[6]);
  double Lmax       = atof(argv[7]);
  int T             = atoi(argv[8]);
  double kbt        = atof(argv[9]);
  double R          = atof(argv[10]);
  sprintf(filename,"%s",argv[11]);
  int posOut        = atoi(argv[12]);
  int forceOut      = atoi(argv[13]);
  int seed          = atoi(argv[14]);
  
  //inits
  double qL = -0.5*N*ci_charge;

  //init arrays
  double *q; 
  double **x; double **xn;
  q = malloc((N+2*Nl)*sizeof(double));
  x = malloc((N+2*Nl)*sizeof(double *));
  xn = malloc((N+2*Nl)*sizeof(double *));

  for(int i=0;i<N+2*Nl;i++)
    {
      x[i] = malloc(3*sizeof(double));
      xn[i] = malloc(3*sizeof(double));
      if(x[i] == NULL || xn[i] == NULL)
	{ printf("Out of memory\n"); }
    }

  //set random seed
  settable(1234567890987654321ULL, seed*1000000000000ULL + 123456123456123456ULL, 362436362436362436ULL, 1066149217761810ULL);


  //open required output files
  char fname [32];
  FILE *pos;
  if (posOut != 0) {
    sprintf(fname, "%s.xyz",filename); pos=fopen(fname,"w");
  }
  char fnameforce [32];
  FILE *force;
  if (forceOut != 0) {
    sprintf(fnameforce, "%s.f",filename); force=fopen(fnameforce,"w");
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
      q[i] = qL/Nl;
      x[i][0] = xn[i][0] = -0.5*R;
      x[i][1] = xn[i][1] = (i-N)*Ly/Nl;
      x[i][2] = xn[i][2] = 0.0;
    }
    else {
      q[i] = qL/Nl;
      x[i][0] = xn[i][0] = 0.5*R;
      x[i][1] = xn[i][1] =  (i-N-Nl)*Ly/Nl;
      x[i][2] = xn[i][2] = 0.0;
    }
  }

  //MC loop
  for(s=0; s < T; s++) {
    ranN = ran_particle(N+2*Nl);
    dx = ran_du(); dy = ran_du(); dz = ran_du();
    step = on_chain(ranN, N) ? 0.2 : 0.6;
    xn[ranN][0] += step*dx;
    xn[ranN][1] += is_endpoint(ranN,N,Nl) ? 0.0 : step*dy;
    xn[ranN][2] += step*dz;
    
    De = delta_u(x,xn,q,ranN,ep,h,htheta,Lmax,N,Nl);
    
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
      fprintf(pos,"N%d Nl%d q%f ep%f h%f Lmax%f T%d kbt%f R%f pos%d\n",N,Nl,ci_charge,ep,h,Lmax,T,kbt,R,posOut);
      for(i=0;i<N+2*Nl;i++) {
        fprintf(pos,"%d %f %f %f\n",on_chain(i,N)?10:1,x[i][0],x[i][1],x[i][2]);
      }
    }

    if(s >= 100000 && forceOut != 0 && s % forceOut == 0) {
      fLx = 0.0;
      fRx = 0.0;
      for(int i=0; i<N+2*Nl; i++) {
        for(int j=0; j< N+2*Nl;j++) {
          if(i==j) { continue; }
          
          if(i >= N && i < N+Nl) {//if i is on left chain, add to fLx
            fLx += lekner_fx(q[i],x[i][0],x[i][1],x[i][2],q[j],x[j][0],x[j][1],x[j][2]);
                //+  rep_fx(ep, x[i][0],x[i][1],x[i][2],x[j][0],x[j][1],x[j][2]);
            //fLx += linked(i,j,N,Nl) ? spring_fx(h, Lmax, x[i][0],x[i][1],x[i][2],x[j][0],x[j][1],x[j][2]) : 0.0;
          }
          else if(i >= N+Nl) {//if i is on the right chain, add to fRx
            fRx += lekner_fx(q[i],x[i][0],x[i][1],x[i][2],q[j],x[j][0],x[j][1],x[j][2]);
                //+  rep_fx(ep, x[i][0],x[i][1],x[i][2],x[j][0],x[j][1],x[j][2]);
            //fRx += linked(i,j,N,Nl) ? spring_fx(h, Lmax, x[i][0],x[i][1],x[i][2],x[j][0],x[j][1],x[j][2]) : 0.0;
          }
        }
      }
      fLavg += fLx*forceOut/(T-100000);
      fRavg += fRx*forceOut/(T-100000);
      fprintf(force, "%f\t%f\n", fLavg, fRavg);
    }
  }

  //close files
  if(posOut != 0) { fclose(pos); }
  if(forceOut != 0) { fclose(force); }

  for(int i=0; i< N+2*Nl; i++)
    {
      free(x[i]);
      free(xn[i]);
    }
  free(x); free(xn); free(q);

  tock(&tt);
}

#include "functions.c"
#define nequil 12000

int main(int argc, char **argv) {
  pca_time tt;
  tick(&tt);
  //gsl_set_error_handler_off();
  
  //variable init
  int s,n;
  double dx, dy, dz, fLx, fRx, kbt, acc_ratio = 0.0;
  double fLprev=0.0, fRprev=0.0;
  int accepted = 0;
  double De;
  char filename[32];
  int i; int ranN;
  bool loop = true;

  //default parameters
  int N             = 16;
  int Nl            = 32;
  double ci_charge  = 1.0;
  double ep         = 1.0;
  double h          = 1.0;
  double htheta     = 1.0;
  double Lmax       = 0.8;
  int T             = 200000;
  double kf         = 1.0;
  double R          = 2.0;
  int posOut        = 1000;
  int forceOut      = 100;
  int seed          = 1;
  double step       = 0.03;

  //import simulation parameters
  for(int i = 1; i < argc;i++) {
    if ( strcmp(argv[i], "-N") == 0 )
      N  = atoi(argv[++i]);
    else if ( strcmp(argv[i], "-Nl") == 0)
      Nl = atoi(argv[++i]);
    else if ( strcmp(argv[i], "-q") == 0 )
      ci_charge  = atof(argv[++i]);
    else if ( strcmp(argv[i], "-e") == 0 )
      ep = atof(argv[++i]);
    else if ( strcmp(argv[i], "-h") == 0 )
      h  = atof(argv[++i]);
    else if ( strcmp(argv[i], "-hth") == 0 )
      htheta = atof(argv[++i]);
    else if ( strcmp(argv[i], "-L") == 0 )
      Lmax   = atof(argv[++i]);
    else if ( strcmp(argv[i], "-T") == 0 )
      T      = atoi(argv[++i]);
    else if ( strcmp(argv[i], "-kf") == 0 )
      kf     = atof(argv[++i]);
    else if ( strcmp(argv[i], "-R") == 0 )
      R      = atof(argv[++i]);
    else if ( strcmp(argv[i], "-po") == 0 )
      posOut = atoi(argv[++i]);
    else if ( strcmp(argv[i], "-fo") == 0 )
      forceOut = atoi(argv[++i]);
    else if ( strcmp(argv[i], "-s") == 0 )
      seed     = atoi(argv[++i]);
    else if ( strcmp(argv[i], "-D") == 0 )
      step     = atof(argv[++i]);
    else if ( strcmp(argv[i], "-f") == 0 )
      strcpy(filename,argv[++i]); 
  }
 
  //inits
  double qL = -0.5*N*ci_charge/(1.0*Nl);
  step       = 0.00075+0.0084*exp((kf-0.02842)/0.0326); //doesn't allow setting a diff stepsize

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
  settable(1234567890987654321ULL, seed*1000 + 123456123456123456ULL, 362436362436362436ULL, 1066149217761810ULL);


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
      q[i]    = ci_charge*0.1;
      x[i][0] = xn[i][0] = ran_xz(R);
      x[i][1] = xn[i][1] = ran_y(Ly);
      x[i][2] = xn[i][2] = ran_xz(0.5*R);
    }
    else if (i<N+Nl) {
      q[i]    = qL*0.1;
      x[i][0] = xn[i][0] = -0.5*R*cos(twoPI*3.0*(i-N)/(Nl-1));
      x[i][1] = xn[i][1] = (i-N)*Ly/Nl;
      x[i][2] = xn[i][2] = -0.5*R*sin(twoPI*3.0*(i-N)/(Nl-1));
    }
    else {
      q[i]    = qL*0.1;
      x[i][0] = xn[i][0] = 0.5*R*cos(twoPI*3.0*(i-N-Nl)/(Nl-1));
      x[i][1] = xn[i][1] =  (i-N-Nl)*Ly/Nl;
      x[i][2] = xn[i][2] = 0.5*R*sin(twoPI*3.0*(i-N-Nl)/(Nl-1));
    }
  }

  //MC loop
  for(s = 0; s < T; s++) {
    for(n = 0; n < N+2*Nl; n++) {
      kbt = kf;
      ranN = ran_particle(N+2*Nl);
      dx = ran_du(); dy = ran_du(); dz = ran_du();

      xn[ranN][0] += step*dx;
      xn[ranN][1] += is_endpoint(ranN,N,Nl) ? 0.0 : step*dy;
      xn[ranN][2] += step*dz;

      if(xn[ranN][1] < 0.0) xn[ranN][1] += Ly;
      else if(xn[ranN][1] > Ly) xn[ranN][1] -= Ly;
      
      De = delta_u(x,xn,q,ranN,ep,h,htheta,Lmax,N,Nl);
      
      if(De < 0 || exp(-De/kbt) > ran_u()) {
        if(s >= nequil) accepted += 1;
        x[ranN][0] = xn[ranN][0];
        x[ranN][1] = xn[ranN][1];
        x[ranN][2] = xn[ranN][2];
      }
      else {
        xn[ranN][0] = x[ranN][0];
        xn[ranN][1] = x[ranN][1];
        xn[ranN][2] = x[ranN][2];
      }
    } //end sweep

    if(s < nequil && s % 100 == 0) {
      for(i=0; i<N+2*Nl;i++) {
        if(i<N) {
          q[i] = ci_charge*(0.1+0.9*s/(1.0*(nequil-100)));
        }
        else {
          q[i] = qL*(0.1+0.9*s/(1.0*(nequil-100)));
        }
      }
#ifdef DEBUG
      printf("qci: %f, qL: %f\n",q[0],q[N]);
#endif
    }
    
    if(posOut != 0 && s >= nequil && s % posOut == 0) {
      fprintf(pos,"%d\n",N+2*Nl);
      fprintf(pos,"N%d Nl%d q%f ep%f h%f Lmax%f T%d kbt%f R%f pos%d\n",N,Nl,ci_charge,ep,h,Lmax,T,kbt,R,posOut);
      for(i=0;i<N+2*Nl;i++) {
        fprintf(pos,"%d %f %f %f\n",on_chain(i,N)?10:1,x[i][0],x[i][1],x[i][2]);
      }
    }

    if(forceOut != 0 && s >= nequil && s % forceOut == 0) {
      fLx = 0.0;
      fRx = 0.0;
      for(int i=0; i<N+2*Nl; i++) {
        for(int j=0; j< N+2*Nl;j++) {
          if(i==j) { continue; }
          
          if(i >= N && i < N+Nl) {//if i is on left chain, add to fLx
            fLx += lekner_fx(q[i],x[i][0],x[i][1],x[i][2],q[j],x[j][0],x[j][1],x[j][2]);
                +  rep_fx(ep, x[i][0],x[i][1],x[i][2],x[j][0],x[j][1],x[j][2]);
            fLx += linked(i,j,N,Nl) ? spring_fx(h, Lmax, x[i][0],x[i][1],x[i][2],x[j][0],x[j][1],x[j][2]) : 0.0;
          }
          else if(i >= N+Nl) {//if i is on the right chain, add to fRx
            fRx += lekner_fx(q[i],x[i][0],x[i][1],x[i][2],q[j],x[j][0],x[j][1],x[j][2]);
                +  rep_fx(ep, x[i][0],x[i][1],x[i][2],x[j][0],x[j][1],x[j][2]);
            fRx += linked(i,j,N,Nl) ? spring_fx(h, Lmax, x[i][0],x[i][1],x[i][2],x[j][0],x[j][1],x[j][2]) : 0.0;
          }
        }
      }
 
      fLprev = fLx;
      fRprev = fRx;

      fprintf(force, "%f\t%f\t%f\n", R, fLx, fRx);
    } //end forceOut

  } //end MC loops

  for(int i=0; i<Nl;i++){
    printf("%f\n",acos((x[N+i][0]-x[N+Nl+i][0])/sqrt((x[N+i][0]-x[N+Nl+i][0])*(x[N+i][0]-x[N+Nl+i][0])+(x[N+i][2]-x[N+Nl+i][2])*(x[N+i][2]-x[N+Nl+i][2]))));
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

  printf("%d\t%d\t\n", accepted, N*(s-nequil));
  tock(&tt);
}

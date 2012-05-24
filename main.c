#include "functions.c"
#define nequil 36000

int main(int argc, char **argv) {
  pca_time tt;
  tick(&tt);
  //gsl_set_error_handler_off();
  
  //variable init
  int s,n;
  double dx, dy, dz, fLx, fRx, kbt, acc_ratio = 0.0;
  double fLprev=0.0, fRprev=0.0;
  double sumsq, rms, xL, zL;
  int accepted = 0;
  double De;
  char filename[64];
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
  int dataOut      = 1000;
  int seed          = 1;
  double step       = 0.03;
  double turns      = 0.0;
  double amp        = 0.0;
  double wv         = 0.0;
  double Ly         = 24.0;
  double ic1[64]    = {0.2813, 0.2401, 0.1369, 0.0202, -0.0601, -0.0799, -0.0498, -0.0044, 0.0221, 0.0157, -0.0124, -0.0381, -0.0432, -0.0268, -0.0040, 0.0074, 0.0, -0.0187, -0.0335, -0.0332, -0.0193, -0.0034, 0.0020, -0.0067, -0.0221, -0.0317, -0.0282, -0.0148, -0.0024, -0.0008, -0.0111, -0.0249, -0.0313, -0.0249, -0.0111, -0.0008, -0.0024, -0.0148, -0.0282, -0.0317, -0.0221, -0.0067, 0.0020, -0.0034, -0.0193, -0.0332, -0.0335, -0.0187, 0.0, 0.0074, -0.0040, -0.0268, -0.0432, -0.0381, -0.0124, 0.0157, 0.0221, -0.0044, -0.0498, -0.0799, -0.0601, 0.0202, 0.1369, 0.2401};
  double ic2[64]    = {-0.0313, -0.0249, -0.0111, -0.0008, -0.0024, -0.0148, -0.0282, -0.0317, -0.0221, -0.0067, 0.0020, -0.0034, -0.0193, -0.0332, -0.0335, -0.0187, 0.0, 0.0074, -0.0040, -0.0268, -0.0432, -0.0381, -0.0124, 0.0157, 0.0221, -0.0044, -0.0498, -0.0799, -0.0601, 0.0202, 0.1369, 0.2401, 0.2813, 0.2401, 0.1369, 0.0202, -0.0601, -0.0799, -0.0498, -0.0044, 0.0221, 0.0157, -0.0124, -0.0381, -0.0432, -0.0268, -0.0040, 0.0074, 0.0, -0.0187, -0.0335, -0.0332, -0.0193, -0.0034, 0.0020, -0.0067, -0.0221, -0.0317, -0.0282, -0.0148, -0.0024, -0.0008, -0.0111, -0.0249 };


 

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
      dataOut = atoi(argv[++i]);
    else if ( strcmp(argv[i], "-s") == 0 )
      seed     = atoi(argv[++i]);
    else if ( strcmp(argv[i], "-D") == 0 )
      step     = atof(argv[++i]);
    else if ( strcmp(argv[i], "-tr") == 0 )
      step     = atof(argv[++i]);
    else if ( strcmp(argv[i], "-f") == 0 )
      strcpy(filename,argv[++i]); 
    else if ( strcmp(argv[i], "-A") == 0 )
      amp      = atof(argv[++i]);
    else if ( strcmp(argv[i], "-wv") == 0 )
      wv       = atof(argv[++i]);
    else if ( strcmp(argv[i], "-Ly") == 0 )
      Ly       = atoi(argv[++i]);
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
  char fname [64];
  FILE *pos;
  if (posOut != 0) {
    sprintf(fname, "%s.xyz",filename); pos=fopen(fname,"w");
  }
  char fnamedata [64];
  FILE *data;
  if (dataOut != 0) {
    sprintf(fnamedata, "%s.dat",filename); data=fopen(fnamedata,"w");
    //print headers
    fprintf(data,"%%seed\tt\tN\tNl\tci_charge\tep\th\ththeta\tLmax\tkf\tR\tturns\tfLx\tfRx\tRl\tA\twv\tA0\n");
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
      x[i][0] = xn[i][0] = -0.5*R*cos(twoPI*turns*(i-N)/(Nl-1))+amp*sin(twoPI*(i-N)*Ly/Nl/wv);
      x[i][1] = xn[i][1] = (i-N)*Ly/Nl;
      x[i][2] = xn[i][2] = -0.5*R*sin(twoPI*turns*(i-N)/(Nl-1));
    }
    else {
      q[i]    = qL*0.1;
      x[i][0] = xn[i][0] = 0.5*R*cos(twoPI*turns*(i-N-Nl)/(Nl-1));
      x[i][1] = xn[i][1] =  (i-N-Nl)*Ly/Nl;
      x[i][2] = xn[i][2] = 0.5*R*sin(twoPI*turns*(i-N-Nl)/(Nl-1));
    }
  }

  //MC loop
  for(s = 0; s < T; s++) {
    for(n = 0; n < N+2*Nl; n++) {
      kbt = kf;
      ranN = ran_particle(N+2*Nl);
      dx = ran_du(); dy = ran_du(); dz = ran_du();

      if( ranN >= N && s < nequil ) { //freeze chain particles before neq
        dx = 0.0; dy = 0.0; dz = 0.0;
      }

      xn[ranN][0] += ranN < N               ? step*dx : step*dx*0.01;
      xn[ranN][1] += is_endpoint(ranN,N,Nl) ? 0.0     : step*dy;
      xn[ranN][2] += ranN < N               ? step*dz : step*dz*0.01;

      if(xn[ranN][1] < 0.0) xn[ranN][1] += Ly;
      else if(xn[ranN][1] > Ly) xn[ranN][1] -= Ly;
      
      De = delta_u(x,xn,q,ranN,ep,h,htheta,Lmax,N,Nl,Ly);
      
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

    if(dataOut != 0 && s >= nequil && s % dataOut == 0) {
      //initialize output vars
      fLx = 0.0;
      fRx = 0.0;
      xL  = 0.0;
      zL  = 0.0;
      rms  = 0.0;
      sumsq= 0.0;

      //calculate forces
      for(int i=0; i<N+2*Nl; i++) {
        for(int j=0; j< N+2*Nl;j++) {
          if(i==j) { continue; }
          
          if(i >= N && i < N+Nl) {//if i is on left chain, add to fLx
            fLx += lekner_fx(q[i],x[i][0],x[i][1],x[i][2],q[j],x[j][0],x[j][1],x[j][2],Ly);
                +  rep_fx(ep, x[i][0],x[i][1],x[i][2],x[j][0],x[j][1],x[j][2],Ly);
            fLx += linked(i,j,N,Nl) ? spring_fx(h, Lmax, x[i][0],x[i][1],x[i][2],x[j][0],x[j][1],x[j][2],Ly) : 0.0;
          }
          else if(i >= N+Nl) {//if i is on the right chain, add to fRx
            fRx += lekner_fx(q[i],x[i][0],x[i][1],x[i][2],q[j],x[j][0],x[j][1],x[j][2],Ly);
                +  rep_fx(ep, x[i][0],x[i][1],x[i][2],x[j][0],x[j][1],x[j][2],Ly);
            fRx += linked(i,j,N,Nl) ? spring_fx(h, Lmax, x[i][0],x[i][1],x[i][2],x[j][0],x[j][1],x[j][2],Ly) : 0.0;
          }
        }
      }
 
      fLprev = fLx;
      fRprev = fRx;

      
      //calculate average position of left line
      for(int i=N; i<N+Nl; i++) {
        xL += x[i][0]/Nl;
        zL += x[i][2]/Nl;
      }

      //calculate rms amplitude of left line
      for(int i=N; i<N+Nl; i++) {
        sumsq += (x[i][0]-xL)*(x[i][0]-xL);
      }
      rms = sqrt(sumsq/Nl);

      fprintf(data,"%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", seed, s, N, Nl, ci_charge, ep, h, htheta, Lmax, kf, R, turns, fLx, fRx, xL, rms, wv,amp);
    } //end dataOut

  } //end MC loops

  /*for(int i=0; i<Nl;i++){
    printf("%f\n",acos((x[N+i][0]-x[N+Nl+i][0])/sqrt((x[N+i][0]-x[N+Nl+i][0])*(x[N+i][0]-x[N+Nl+i][0])+(x[N+i][2]-x[N+Nl+i][2])*(x[N+i][2]-x[N+Nl+i][2]))));
  }*/

  //close files
  if(posOut != 0) { fclose(pos); }
  if(dataOut != 0) { fclose(data); }

  for(int i=0; i< N+2*Nl; i++)
    {
      free(x[i]);
      free(xn[i]);
    }
  free(x); free(xn); free(q);

  printf("%d\t%d\t\n", accepted, N*(s-nequil));
  tock(&tt);
}

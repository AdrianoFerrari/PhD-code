#include "tests.c"

int main(int argc, char **argv) {
  //run_tests(); 

  //import simulation parameters
  int N = atoi(argv[1]);
  int s; 
  double maxR = 16.0; double maxY = 24.0;
  double ci_charge = atof(argv[2]);
  double kbt;
  double ep = atof(argv[3]);
  double qL = -0.5*N*ci_charge;
  char filename[32];
  int T = atoi(argv[4]);
  double kbt0 = atof(argv[5]);
  double kbtf = atof(argv[6]);
  double R = atof(argv[7]);
  double A = atof(argv[8]);
  double lambda = atof(argv[9]);
  sprintf(filename,"%s",argv[10]);
  int growthOut = atoi(argv[11]);
  int forceOut = atoi(argv[12]);
  int posOut = atoi(argv[13]);

  //declare variables other variables
  int i;  double De; int ranN;
  double q[N];
  double x[N][3];
  double xn[N][3];
  
  //initialize MPI
  int node; int size; ULL idum;
  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&node);
  idum = 1234567890987654321ULL + node*100000000000000ULL;
  settable(idum, 123456123456123456ULL, 362436362436362436ULL, 1066149217761810ULL);

  //open required output files
  char fname [32];
  FILE *fp; FILE *force; FILE *pos;
  if (growthOut != 0) { sprintf(fname, "%s.%d.growth",filename,node); fp=fopen(fname,"w"); }
  if (forceOut != 0) { sprintf(fname, "%s.%d.force",filename,node); force=fopen(fname,"w"); }
  if (posOut != 0) { sprintf(fname, "%s.%d.xyz",filename,node); pos=fopen(fname,"w"); }

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
  for(s=0; s < T/size; s++) {
    ranN = ran_particle(N);
    xn[ranN][0] = ran_xz(maxR);
    xn[ranN][1] = ran_y(maxY);
    xn[ranN][2] = ran_xz(maxR);
    
    kbt = (kbtf - kbt0)*s/(T/size-1) + kbt0;
    
    if(go_ahead(energy_difference(x,xn,q,N,ranN,ep,qL,-0.5*R,0.0,A,0.5*R,0.0,A,lambda), kbtf)) {
      x[ranN][0] = xn[ranN][0];
      x[ranN][1] = xn[ranN][1];
      x[ranN][2] = xn[ranN][2];
    }
    else {
      xn[ranN][0] = x[ranN][0];
      xn[ranN][1] = x[ranN][1];
      xn[ranN][2] = x[ranN][2];
    }
    
    if(growthOut != 0 && s % growthOut == 0) {
      fprintf(fp,"%f\n",growth_rate(q,x,N,qL,-0.5*R,0.0,0.5*R,0.0,ep,A,lambda,23,7));
    }
    if(forceOut != 0 && s % forceOut == 0) {
      fprintf(force,"%f\n", xforce(q,x,N,qL,-0.5*R,0.0,0.5*R,0.0,ep,A,lambda,23,7));
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
  if(growthOut != 0) { fclose(fp); }
  if(forceOut != 0) { fclose(force); }
  if(posOut != 0) { fclose(pos); }
  
  //close MPI
  MPI_Finalize();
}

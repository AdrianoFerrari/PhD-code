#include "tests.c"

int main(int argc, char **argv) {
  //run_tests(); 

  //import simulation parameters
  int N = 64;
  int s; 
  double maxR = 16.0; double maxY = 24.0;
  double ci_charge = 2.0;
  double kbt;
  double ep = 1.0;
  double qL = -5.333333;
  char filename[32];
  int T = atoi(argv[1]);
  double kbt0 = atof(argv[2]);
  double kbtf = atof(argv[3]);
  double R = atof(argv[4]);
  double A = atof(argv[5]);
  double lambda = atof(argv[6]);
  sprintf(filename,"%s",argv[7]);
  int outputrate = atoi(argv[8]);

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
  FILE *fp;
  char fullfilename [32];
  sprintf(fullfilename, "%s.%d",filename,node);
  fp=fopen(fullfilename,"w");
  FILE *pos;
  char posfilename [32];
  sprintf(posfilename, "%s.%d.xyz",filename,node);
  pos=fopen(posfilename,"w");

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
    ranN = (int)floor(rand()/RAND_MAX);
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
    
    if(s % outputrate == 0) {
      fprintf(fp,"%f, %f\n",growth_rate(q,x,N,qL,-0.5*R,0.0,0.5*R,0.0,ep,A,lambda,23,7),kbt);
      
      fprintf(pos,"%d\n",N);
      fprintf(pos,"somestuff\n");
      for(i=0;i<N;i++) {
	fprintf(pos,"0 %f %f %f\n",x[i][0],x[i][1],x[i][2]);
      }
    }
  }

  //close files
  fclose(fp);
  
  //close MPI
  MPI_Finalize();
}

#include "tests.c"

int main() {
  //run_tests(); 

  //import simulation parameters
  int N = 16;
  double maxR = 16.0; double maxY = 24.0;
  double ci_charge = 2.0;
  double kbT = 1.0;
  double ep = 1.0;
  double qL = -10; double R = 1.2; double A = 0.1; double lambda = 6.0;
  double De;
  //declare variables
  int i;
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
  fp=fopen("testout1","w");

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
  for(s=0; s < T/size;s++) {
    ranN =(int)floor(N*UNI);
    xn[ranN][0] = ran_xz(maxR);
    xn[ranN][1] = ran_y(maxY);
    xn[ranN][2] = ran_xz(maxR);
    
    if(go_ahead(energy_difference(x,xn,q,N,ranN,ep,qL,-0.5*R,0.0,A,0.5*R,0.0,A,lambda))) {
      x[ranN][0] = xn[ranN][0];
      x[ranN][1] = xn[ranN][1];
      x[ranN][2] = xn[ranN][2];
    }
    else {
      xn[ranN][0] = x[ranN][0];
      xn[ranN][1] = x[ranN][1];
      xn[ranN][2] = x[ranN][2];
    }
    
    if(s % 10000 == 0) {
      fprintf(fp,"%d",growth_rate(q,x,N,qL,-0.5*R,0.0,0.5*R,0.0,ep,A,lambda,23,7));
    }
  }

  //close files
  fclose(fp);
  
  //close MPI
  MPI_Reduce(&myenergy,&energy,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Finalize();
}

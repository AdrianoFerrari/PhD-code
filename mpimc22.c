#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "bessels.h"

#define PI 3.14159265358979323846
#define twoPI 6.283185307179586
#define znew (z=36969*(z&65535)+(z>>16))
#define wnew (w=18000*(w&65535)+(w>>16))
#define MWC ((znew<<16)+wnew )
#define SHR3 (jsr^=(jsr<<17), jsr^=(jsr>>13), jsr^=(jsr<<5))
#define CONG (jcong=69069*jcong+1234567)
#define FIB ((b=a+b),(a=b-a))
#define KISS ((MWC^CONG)+SHR3)
#define LFIB4 (c++,t[c]=t[c]+t[UC(c+58)]+t[UC(c+119)]+t[UC(c+178)])
#define SWB (c++,bro=(x<y),t[c]=(x=t[UC(c+34)])-(y=t[UC(c+19)]+bro))
#define UNI (KISS*5.421010862e-20)
#define VNI ((long) KISS)*1.084202172e-19
#define UC (unsigned char) /*a cast operation*/
typedef unsigned long UL;
/* Global static variables: */
static UL z=362436069, w=521288629, jsr=123456789, jcong=380116160;
static UL a=224466889, b=7584631, t[256];
/* Use random seeds to reset z,w,jsr,jcong,a,b, and the table
t[256]*/
static UL x=0,y=0,bro; static unsigned char c=0;

typedef struct PARAMS{
unsigned int T;
double kbTF;
double amp;
double lambda;
double R;
char base[32];
double kbT0;
int seed;
int coords_out;
int force_out;
int growth_out;
double sigma;
double rljEps;
int Ns;
} PARAMS;

const int Nx = 64;//lattice points
const int Ny = 96;//lattice points
const int Nz = 10;//lattice points
const double Lx = 16.0;
const double Ly = 24.0;
const double uy = 0.041666666666;
const int N = 64;
//const double sigma = 0.45;
//const double rljEps = 1.0;
const double linCharge = -5.33333333;
const int M = 7;//Lekner terms
const double besskx[32]={0.0000161803, 0.0000261803, 0.0000423607,0.000068541, 0.000110902,0.000179443, 0.000290344, 0.000469787, 0.000760132, 0.00122992,0.00199005, 0.00321997, 0.00521002, 0.00842999, 0.01364, 0.02207,0.03571, 0.05778, 0.09349,0.15127, 0.24476, 0.39603, 0.64079,1.03682, 1.67761, 2.71443, 4.39204, 7.10647, 11.4985, 18.605, 30.1035,48.7085};
const double bessk0y[32]={11.1476,10.6664,10.1852,9.70401,9.2228,8.74159,8.26037,7.77916,7.29795,6.81674,6.33553,5.85433,5.37315,4.892,4.41093,3.93007,3.44967,2.97036,2.49345,2.02184,1.56137,1.12325,0.726727,0.399538,0.170257,0.0484295,0.00721246,0.00037911,0.0000037108,0.00000000240083,0.0000000000000191954, 0.0};
const double bessk1y[32]={61803.4,38196.6,23606.8,14589.8,9016.99,5572.81,3444.18,2128.62,1315.56,813.057,502.493,310.552,191.923,118.601,73.2803,45.2615,27.9328,17.2068,10.5566,6.42055,3.83546,2.21072,1.19023,0.5657,0.21594,0.0567286,0.00799482,0.000404955,0.00000386892,0.00000000246453,0.0000000000000195117,0.0};

static double CalculateE(double q[],double k[][],PARAMS pars);
static double CalculateEn(double q[],double k[][],PARAMS pars, int ranN);
static double CalculateE_n(double q[],double k[][],PARAMS pars,int ranN);

static double LeknerPotentialE(double q[],double k[][]);
static double LeknerPotentialEn(double q[],double k[][],int ranN);
static double LeknerPotentialE_n(double q[],double k[][],int ranN);
static double Lekner(double qi,double qj, double xi,double yi,double zi,double xj,double yj,double zj);

static double RLJPotentialE(double k[][],PARAMS pars);
static double RLJPotentialEn(double k[][],int ranN,PARAMS pars);
static double RLJPotentialE_n(double k[][],int ranN,PARAMS pars);
static double RLJ( double xi,double yi,double zi,double xj,double yj,double zj,PARAMS pars);

static double LinePotentialE(double q[],double k[][],PARAMS pars, double x0,double z0);
static double LinePotentialEn(double q[],double k[][],PARAMS pars, double x0,double z0,int ranN);
static double LinePotentialE_n(double q[],double k[][],PARAMS pars, double x0,double z0,int ranN);
static double Line(double qi,double xi,double yi,double zi,PARAMS pars,double x0,double z0);

static double LineRLJPotentialE(double k[][],double x0,double z0,PARAMS pars);
static double LineRLJPotentialEn(double k[][],double x0,double z0,PARAMS pars,int ranN);
static double LineRLJPotentialE_n(double k[][],double x0,double z0, PARAMS pars,int ranN);
static double LineRLJ(double xi,double yi,double zi,PARAMS pars,double x0,double z0);

static double ForceCoul_LineParticle(float q, double x, double y, double z, double x0, double z0);
static double ForceRLJ_LineParticle(double x, double y, double z, double x0, double z0,PARAMS pars);
static double ForceLine_Line(double x0, double z0, double x1, double z1);
static double GrowthRate(double s, float q, double x, double y, double z, double x0, double z0, PARAMS pars);
static double AverageGrowthRate(double q[],double k[][],PARAMS pars);

void settable(UL i1,UL i2,UL i3,UL i4,UL i5, UL i6)
{ int i; z=i1;w=i2,jsr=i3; jcong=i4; a=i5; b=i6;
for(i=0;i<256;i=i+1) t[i]=KISS;
}

int main(int argc,char **argv)
{
 //initialize variables
   int node;
   int size;
   unsigned long idum;
   struct PARAMS pars;
   pars.T=atoi(argv[1]);
   pars.kbTF=atof(argv[2]);
   pars.amp=atof(argv[3]);
   pars.lambda=atof(argv[4]);
   pars.R=atof(argv[5]);
   sprintf(pars.base,"%s",argv[6]);
   pars.kbT0=atof(argv[7]);
   double kbT;
   double kbTlist[10]= {5.0,4.0,3.0,2.0,1.5,1.0,0.5,0.3,0.1,0.05};
   double q[N];
   double k[3][N];
   double kTemp[3][N];
   int ranN;
   double ranR;
   double ranPhi;
   double ranY;
   double De;
   double myenergy=0.0;
   double energy=0.0;
   double dx=1;
   double accepted=0.0;
   double accValsTested=0.0;
   pars.seed = atoi(argv[8]); 
   pars.coords_out = atoi(argv[9]);
   pars.force_out = atoi(argv[10]);
   pars.growth_out = atoi(argv[11]);
   pars.sigma = atof(argv[12]);
   pars.rljEps = atof(argv[13]);
   pars.Ns = atoi(argv[14]);

//set charges
   for(int i=0; i<N;i++)
   {q[i]=2.0;}

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&size);
   MPI_Comm_rank(MPI_COMM_WORLD,&node);
   
   //Each node has diff initial conditions
   idum=pars.seed+node*100000000;
   
   //const char base[]="cRun.";
   char filename [32];
   sprintf(filename, "%s%d",pars.base, node);

   FILE *fp;   
   if(pars.coords_out!=0){
   fp=fopen(filename,"w");
   }

   char forceFilename [32];
   sprintf(forceFilename, "for_%s%d",pars.base, node);

   FILE *ffp;
   if(pars.force_out!=0){
   ffp=fopen(forceFilename,"w");
   }

   char grwthfilename [32];
   sprintf(grwthfilename, "grw_%s%d",pars.base, node);
   FILE *grwth;
   if(pars.growth_out!=0){
   grwth=fopen(grwthfilename,"w");
   }

   settable(362436069,521288629,idum,380116160,224466889,7584631);

   for(int i=0;i<N;i++)
   {
      ranR=sqrt(Lx*Lx*UNI);
      ranPhi=2*PI*UNI;
      ranY=Ly*UNI;

      k[0][i]=ranR*cos(ranPhi);
      k[1][i]=ranY;
      k[2][i]=ranR*sin(ranPhi);
   }//NOTE: did not exclude initializing charges within cyl!!

   //Main MC loop
   double E=CalculateE(q,k,pars);//initial energy calculation
   //printf("%f\n",E);   
   int Gdenom = 0;
   double avgG = 0.0;

   for(int s=0;s<pars.T/size;s++)
   {
	ranN=(int)floor(N*UNI);
	for(int i=0;i<N;i++)
	{
        kTemp[0][i]=k[0][i];
        kTemp[1][i]=k[1][i];
        kTemp[2][i]=k[2][i];
	}

        kbT=(pars.kbTF-pars.kbT0)*s/(pars.T-1)+pars.kbT0;        
	
	//test move
        ranR=sqrt(Lx*Lx*UNI);
        ranPhi=2*PI*UNI;
        ranY=Ly*UNI;
        
	
	kTemp[0][ranN]=ranR*cos(ranPhi);
	kTemp[1][ranN]=ranY;
        kTemp[2][ranN]=ranR*sin(ranPhi);

	double En=CalculateEn(q,k,pars,ranN);
	double Enp=CalculateEn(q,kTemp,pars,ranN);
	De=Enp-En;
         
        if(node==0){
accValsTested++;
//printf("%f\n",De);
	}

    if(De<0)
	{
	   k[0][ranN]=kTemp[0][ranN];
	   k[1][ranN]=kTemp[1][ranN];
           k[2][ranN]=kTemp[2][ranN];
	   if(node==0){accepted++;}
	   double E_n=CalculateE_n(q,kTemp,pars,ranN);
	   E = E_n+Enp;
	}
	else if(exp(-De/kbT)>=UNI)
	{
	   k[0][ranN]=kTemp[0][ranN];
           k[1][ranN]=kTemp[1][ranN];
           k[2][ranN]=kTemp[2][ranN];
	   if(node==0){accepted++;}
	   double E_n=CalculateE_n(q,kTemp,pars,ranN);
	   E = E_n+Enp;
	}

	if(pars.coords_out!=0 && s%pars.coords_out==0)//sanity check, output coords
	{
		fprintf(fp,"%d\n",N);
		fprintf(fp,"%d %f %f %f %f %f\n",pars.T,kbT,pars.amp,pars.lambda,pars.R,pars.kbT0);
		for(int i=0;i<N;i++)
		{
		   fprintf(fp,"0 %f %f %f\n",k[0][i],k[1][i],k[2][i]);
		}
	}

	if(pars.force_out!=0 && s%pars.force_out ==0)
	{
	 double FrljL=0.0;
	 double FrljR=0.0;
	 double FplL=0.0;
	 double FplR=0.0;
		
	 for (int j=0;j<N;j++)
	 {
	   FrljL+=ForceRLJ_LineParticle(k[0][j],k[1][j],k[2][j],-0.5*pars.R,0.0,pars);
	   FrljR+=ForceRLJ_LineParticle(k[0][j],k[1][j],k[2][j],0.5*pars.R,0.0,pars);	
	   FplL+=ForceCoul_LineParticle(q[j],k[0][j],k[1][j],k[2][j],-0.5*pars.R,0.0);
	   FplR+=ForceCoul_LineParticle(q[j],k[0][j],k[1][j],k[2][j],0.5*pars.R,0.0);
	 }
	 fprintf(ffp,"%f,%f,%f,%f\n",FplL,FplR,FrljL,FrljR);
	}
	
	if(pars.growth_out!=0 && s%pars.growth_out==0)
	{
         avgG += AverageGrowthRate(q,k,pars);
	 Gdenom++;
	}
   }//end of main s loop
   
   if(pars.growth_out!=0)
   {
   avgG = avgG/(1.0*Gdenom);
   fprintf(grwth,"%f\n",avgG);
   }
 
MPI_Reduce(&myenergy,&energy,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

//Output
if(node==0)
{
   printf("Energy: %f\nAcceptance Rate: %f\n",energy,accepted/accValsTested);   
}           
   if(pars.coords_out!=0){fclose(fp);}
   if(pars.force_out!=0){fclose(ffp);}
   if(pars.growth_out!=0){fclose(grwth);}
   MPI_Finalize();
}//end of main
static double CalculateE(double q[],double k[][],PARAMS p)
{
   return LeknerPotentialE(q,k)+RLJPotentialE(k,p)+LinePotentialE(q,k,p,-0.5*p.R,0)+LinePotentialE(q,k,p,0.5*p.R,0)+LineRLJPotentialE(k,-0.5*p.R,0,p)+LineRLJPotentialE(k,0.5*p.R,0,p);
}
static double CalculateEn(double q[],double k[][],PARAMS p, int ranN)
{
	return LeknerPotentialEn(q,k,ranN)+RLJPotentialEn(k,ranN,p)+LinePotentialEn(q,k,p,-0.5*p.R,0,ranN)+LinePotentialEn(q,k,p,0.5*p.R,0,ranN)+LineRLJPotentialEn(k,-0.5*p.R,0,p,ranN)+LineRLJPotentialEn(k,0.5*p.R,0,p,ranN);
}
static double CalculateE_n(double q[],double k[][],PARAMS p, int ranN)
{
	return LeknerPotentialE_n(q,k,ranN)+RLJPotentialE_n(k,ranN,p)+LinePotentialE_n(q,k,p,-0.5*p.R,0,ranN)+LinePotentialE_n(q,k,p,0.5*p.R,0,ranN)+LineRLJPotentialE_n(k,-0.5*p.R,0,p,ranN)+LineRLJPotentialE_n(k,0.5*p.R,0,p,ranN);
}

static double LeknerPotentialE(double q[],double x[3][N])
{
	double pe=0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (i != j)
			{
				pe += Lekner(q[i],q[j],x[0][i],x[1][i],x[2][i],x[0][j],x[1][j],x[2][j]);
			}
		}
	}
    return pe;
	
}
static double LeknerPotentialEn(double q[],double x[3][N], int ranN)
{
	double pe=0;
	int i=ranN;
		for (int j = 0; j < N; j++)
		{
			if (i != j)
			{
				pe += Lekner(q[i],q[j],x[0][i],x[1][i],x[2][i],x[0][j],x[1][j],x[2][j]);
			}
		}

    return pe;
	}
static double LeknerPotentialE_n(double q[],double x[3][N], int ranN)
{
	double pe=0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (i != j && i!=ranN && j!=ranN)
			{
				pe += Lekner(q[i],q[j],x[0][i],x[1][i],x[2][i],x[0][j],x[1][j],x[2][j]);
			}
		}
	}
    return pe;
	
}
static double Lekner(double qi,double qj, double xi,double yi,double zi,double xj,double yj,double zj)
{
	double rxz = sqrt((xj - xi) * (xj - xi) + (zj - zi) * (zj - zi));
	double y = (yj - yi);
	double vi = 0;
for (int n = 1; n <= M; n++)
	{
	vi += rxz <= 0 ? 0 : 
	4.0*qi*qj*cos(twoPI*n*y*uy)*bessk0(twoPI*n*rxz*uy)*uy;
	}
	return vi - 2.0*qi*qj*log(rxz)*uy;
}


static double RLJPotentialE(double x[3][N],PARAMS p)
{
        double r;
        double pe = 0;			
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (i != j)
                {
                        pe += RLJ(x[0][i],x[1][i],x[2][i],x[0][j],x[1][j],x[2][j],p);
                }
            }
        }
        return pe;
}
static double RLJPotentialEn(double x[3][N],int ranN,PARAMS p)
{
        double r;
        double pe = 0;			
		int i=ranN;
		for (int j = 0; j < N; j++)
		{
			if (i != j)
			{
					pe += RLJ(x[0][i],x[1][i],x[2][i],x[0][j],x[1][j],x[2][j],p);
			}
		}
        return pe;
}
static double RLJPotentialE_n(double x[3][N], int ranN,PARAMS p)
{
        double r;
        double pe = 0;			
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (i != j && i!=ranN && j!=ranN)
                {
                        pe += RLJ(x[0][i],x[1][i],x[2][i],x[0][j],x[1][j],x[2][j],p);
                }
            }
        }
        return pe;
}
static double RLJ( double xi,double yi,double zi,double xj,double yj,double zj,PARAMS p)
{
double r = sqrt((xj - xi) * (xj - xi) +(yj - yi) * (yj - yi)+ (zj - zi) * (zj - zi));
if (r < 1.122 * p.sigma)
	return 4.0 * p.rljEps * (pow(p.sigma / r,12.0) - pow(p.sigma / r,6.0) + 0.25);
else
	return 0.0;
}


static double LinePotentialE(double q[N],double x[3][N],PARAMS p,double x0,double z0)
{
    double pe = 0;
	
    for (int i = 0; i < N; i++)
    {
        pe += Line(q[i],x[0][i],x[1][i],x[2][i],p,x0,z0);
    }
    return pe;
}
static double LinePotentialEn(double q[N],double x[3][N],PARAMS p,double x0,double z0, int ranN)
{
    double pe = 0;
	
    int i=ranN;
    pe += Line(q[i],x[0][i],x[1][i],x[2][i],p,x0,z0);
    return pe;
}
static double LinePotentialE_n(double q[N],double x[3][N],PARAMS p,double x0,double z0, int ranN)
{
    double pe = 0;
	
    for (int i = 0; i < N; i++)
    {
	if(i!=ranN)
        pe += Line(q[i],x[0][i],x[1][i],x[2][i],p,x0,z0);
    }
    return pe;
}
static double Line(double qi,double xi,double yi,double zi,PARAMS p,double x0,double z0)
{
	double rxz = sqrt((xi - x0) * (xi - x0) + (zi - z0) * (zi - z0));
	if(rxz<=0)
		return 0.0;
	else
	return -2.0*qi*linCharge*log(rxz)+
		4.0*p.amp*PI*qi*xi*linCharge*bessk1(twoPI*rxz/p.lambda)*sin(twoPI*yi/p.lambda)/rxz;
}


static double LineRLJPotentialE(double x[3][N],double x0,double z0,PARAMS p)
{
    double r;
    double pe = 0;
	
    for (int i = 0; i < N; i++)
    {
	pe += LineRLJ(x[0][i],x[1][i],x[2][i],p,x0,z0);
    }
    return pe;
}
static double LineRLJPotentialEn(double x[3][N],double x0,double z0,PARAMS p, int ranN)
{
    double r;
    double pe = 0;
	
    int i=ranN;
	pe += LineRLJ(x[0][i],x[1][i],x[2][i],p,x0,z0);
    return pe;
}
static double LineRLJPotentialE_n(double x[3][N],double x0,double z0,PARAMS p, int ranN)
{
    double r;
    double pe = 0;
	
    for (int i = 0; i < N; i++)
    {
	if(i!=ranN)
		pe += LineRLJ(x[0][i],x[1][i],x[2][i],p,x0,z0);
    }
    return pe;
}
static double LineRLJ(double xi,double yi,double zi,PARAMS p,double x0,double z0)
{
double xa = x0+p.amp*sin(twoPI*yi/p.lambda);
double rxz = sqrt((xi-xa)*(xi-xa) + (zi-z0)*(zi-z0));
if (rxz < 1.122 * p.sigma)
	return 4.0*p.rljEps*(pow(p.sigma/rxz,12.0) - pow(p.sigma/rxz,6.0) + 0.25);
else
	return 0;
}

static double ForceCoul_LineParticle(float q, double x, double y, double z, double x0, double z0)
{
	double rxz = sqrt((x - x0) * (x - x0) + (z - z0) * (z - z0));
	if (x0>x)
		return (2.0*q*linCharge/rxz)*sqrt(1-z*z/(rxz*rxz));
	else
		return -1.0*(2.0*q*linCharge/rxz)*sqrt(1-z*z/(rxz*rxz));		
}
static double ForceRLJ_LineParticle(double x, double y, double z, double x0, double z0,PARAMS p)
{
	double rxz = sqrt((x - x0) * (x - x0) + (z - z0) * (z - z0));
	double sigma6=pow(p.sigma,6);
	double rxz6 = pow(rxz,6);
	if(rxz < 1.122* p.sigma)
		return -24.0*(x-x0)*p.rljEps*sigma6*(2.0*sigma6-rxz6)/(rxz6*rxz6*rxz*rxz);
	else
		return 0.0;
}
static double ForceLine_Line(double x0,double z0,double x1,double z1)
{
	double rxz=sqrt((x0-x1)*(x0-x1)+(z0-z1)*(z0-z1));
	if (rxz>0)
		return linCharge*linCharge*Ly/(4.0*rxz);
	else
		return 0.0;
}

static double GrowthRate(double s, float q, double x, double y, double z, double x0, double z0, PARAMS p)
{
	double dx = (x - x0 - p.amp*sin(twoPI*s*Ly/p.lambda));
	double rxz = sqrt(dx*dx  + (z - z0) * (z - z0));
        double G = 0.0;
	double vi = 0;
	for (int n = 1; n <= M; n++)
        {
        	vi += rxz <= 0 ? 0 :	
8.0*PI*linCharge*q*n*cos(twoPI*n*(y-s*Ly)*uy)*bessk1(twoPI*n*rxz*uy)*uy/(1.0*p.Ns);
        }
        G = vi + 2.0*linCharge*q/(rxz*p.Ns);
	return G/(p.amp*sin(twoPI*s*Ly/p.lambda)*sqrt(1+(z-z0)*(z-z0)/(dx*dx)));
}
static double AverageGrowthRate(double q[N],double x[3][N],PARAMS p)
{
	double avgGL = 0;
	double avgGR = 0;
	for (int j = 1; j < p.Ns; j++)
    	{
	  for (int i = 0; i < N; i++){
	     	avgGL += GrowthRate(j/(1.0*p.Ns),q[i],x[0][i],x[1][i],x[2][i],-p.R*0.5,0,p);
          	avgGR += GrowthRate(j/(1.0*p.Ns),q[i],x[0][i],x[1][i],x[2][i],p.R*0.5,0,p);
	  }
	}
	return (avgGL+avgGR)*0.5;
}

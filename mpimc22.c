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
int T;
double kbTF;
double amp;
double lambda;
double R;
char base[32];
double kbT0;
int seed;
int coords_out;
double sigma;
double rljEps;
} PARAMS;

const int Nx = 64;//lattice points
const int Ny = 96;//lattice points
const int Nz = 10;//lattice points
const int Ns = 31;//growth rate samplings
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
static double AverageGrowthRate(double q[],double k[][],PARAMS pars, int Ns);

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
   pars.sigma = atof(argv[10]);
   pars.rljEps = atof(argv[11]);

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
   if(pars.coords_out==1){
   fp=fopen(filename,"w");
   }

   char forceFilename [32];
   sprintf(forceFilename, "for_%s%d",pars.base, node);

   FILE *ffp;   
   ffp=fopen(forceFilename,"w");

   char grwthfilename [32];
   sprintf(grwthfilename, "grw_%s%d",pars.base, node);
   FILE *grwth;
   grwth=fopen(grwthfilename,"w");

   settable(362436069,521288629,idum,380116160,224466889,7584631);
   if(1==1)
   {
     for(int i=0;i<N;i++)
     {
      ranR=sqrt(Lx*Lx*UNI);
      ranPhi=2*PI*UNI;
      ranY=Ly*UNI;

      k[0][i]=ranR*cos(ranPhi);
      k[1][i]=ranY;
      k[2][i]=ranR*sin(ranPhi);
     }//NOTE: did not exclude initializing charges within cyl!!
   }
   else
   {
     double vals[192]={-11.851990,19.018087,-6.088156,14.543560,20.506149,2.224342,-5.677497,19.916444,-13.854407,-12.669699,7.948299,-7.098848,-0.258067,8.424614,1.594937,-10.297615,16.341029,-10.455664,-0.890589,7.698188,-6.442893,-0.655380,19.197009,10.152294,6.340795,4.610594,-2.723341,5.318855,1.827038,-2.495710,-5.415150,15.590445,14.789276,-11.037197,13.220661,10.246705,-5.754359,10.187361,-7.719106,8.704758,19.200300,-5.128541,6.997408,12.904411,0.907965,-12.210664,17.945483,1.090790,-2.498156,0.634381,4.101275,-11.758001,11.111510,-10.202697,8.554182,13.522557,-0.060956,-0.265857,4.520847,-1.137483,-5.813194,7.444360,-1.103674,6.690192,18.438627,0.961493,-14.583496,20.619232,-4.999082,-3.181695,10.898322,-0.311461,3.706256,3.152496,0.002670,-3.466341,19.448986,-0.293293,-3.837955,7.808658,-0.043926,-3.393624,13.908907,-0.209863,3.669055,8.263646,-0.118999,3.667535,21.449135,0.154077,-2.778444,21.379107,-0.181331,-3.168982,18.323821,0.300616,3.682683,19.211363,0.158413,3.669006,23.755456,-0.191563,-3.673541,17.661391,-0.526439,2.761778,0.788354,-1.899462,3.779273,14.160939,-0.222586,3.690215,21.690651,-0.189479,3.432776,3.798038,-0.214940,-3.797590,1.659790,-0.063146,3.112118,7.493183,0.335369,3.309264,14.081664,0.098270,-3.663559,8.429304,0.159685,-3.450197,8.909661,0.188714,-3.575783,23.165283,0.175137,3.596263,11.996036,-0.185750,-3.444374,5.611344,-0.164676,-3.420105,2.750252,0.200047,-3.297188,14.004038,0.177179,3.420026,6.908852,0.226455,-3.186457,15.982058,-0.098544,-3.368929,5.506680,0.152180,3.627858,17.265310,0.144750,-3.667487,10.861763,0.109997,3.312923,19.904080,0.049121,-3.326798,2.672055,0.070727,3.020681,10.439815,0.027040,-3.275966,23.166632,0.028703,-3.533727,16.835648,0.222726,-3.549321,3.974895,0.245135,-3.472028,11.810684,0.238604,3.566448,2.445299,0.308980,3.636652,20.809467,-0.190696,3.725081,5.400863,0.016798};
     for(int i=0;i<N;i++)
     {
         k[0][i]=vals[3*i];
         k[1][i]=vals[3*i+1];
         k[2][i]=vals[3*i+2];
     }
   }

   //Main MC loop
   double E=CalculateE(q,k,pars);//initial energy calculation
   //printf("%f\n",E);
   
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

	if(s%1000==0)//sanity check, output coords
	{
	     if(pars.coords_out==1){
		fprintf(fp,"%d\n",N);
		fprintf(fp,"%d %f %f %f %f %f\n",pars.T,kbT,pars.amp,pars.lambda,pars.R,pars.kbT0);
 		for(int i=0;i<N;i++)
		{
		   fprintf(fp,"0 %f %f %f\n",k[0][i],k[1][i],k[2][i]);
		}
      	     }

		if(s%4000 ==0)
		{
		 //calculate forces and printout
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

	         double g1 = GrowthRate(0.04166666,2.0,k[0][0],k[1][0],k[2][0],-0.5*pars.R,0.0,pars);
                 fprintf(grwth,"%f,%f,%f,%f= %f\n",0.04166666,k[0][0],k[1][0],k[2][0],g1);
		}
	}

   }
 
MPI_Reduce(&myenergy,&energy,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

//Output
if(node==0)
{
   printf("Energy: %f\nAcceptance Rate: %f\n",energy,accepted/accValsTested);   
}           
   if(pars.coords_out==1){
   fclose(fp);}
   fclose(ffp);
   fclose(grwth);
   MPI_Finalize();

}
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
	double rxzS = sqrt((x - x0) * (x - x0) + (z - z0) * (z - z0));
        double G = 0.0;
	double vi = 0;
	for (int n = 1; n <= M; n++)
        {
        	vi += rxz <= 0 ? 0 :	
8.0*PI*linCharge*q*n*cos(twoPI*n*(y-s*Ly)*uy)*(bessk1(twoPI*n*rxz*uy)-bessk1(twoPI*n*rxzS*uy))*uy/(1.0*Ns);
        }
        G = vi + 2.0*linCharge*q*(1/rxz-1/rxzS)/(1.0*Ns);
	return G/(p.amp*sin(twoPI*s*Ly/p.lambda)*sqrt(1+(z-z0)*(z-z0)/(dx*dx)));
}
static double AverageGrowthRate(double q[N],double x[3][N],PARAMS p, int Ns)
{
	double avgGL = 0;
	double avgGR = 0;
	for (int j = 0; j < Ns; j++)
    	{
	  for (int i = 0; i < N; i++){
		avgGL += GrowthRate(j/Ns,q[i],x[0][i],x[1][i],x[2][i],-p.R/2,0,p);
		avgGR += GrowthRate(j/Ns,q[i],x[0][i],x[1][i],x[2][i],p.R/2,0,p);
	  }
	}
	return (avgGL+avgGR)*0.5;
}

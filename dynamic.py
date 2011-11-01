import sys, random, time
from scipy import sqrt, cos, pi, log, exp
from scipy.special import k0
from numpy import zeros, array
from copy import deepcopy
from mpi4py import MPI

maxR = 16.0
Ly = 24.0
uy = 1.0/Ly
M = 7

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def dist(y1,y0):
	if abs(y1-y0) > Ly/2:
		return Ly/2-abs(y1-y0)
	else:
		return abs(y1-y0)

def lekner_u(q0,x0,y0,z0,q1,x1,y1,z1):
	rxz = sqrt((x1-x0)**2 + (z1-z0)**2)
	y 	= dist(y1,y0)
	vi	= 0.0
	for n in range(1,M+1):
		vi += 4.0*q0*q1*cos(2*pi*n*y*uy)*k0(2*pi*n*rxz*uy)*uy if not rxz == 0 else 0
	return vi - 2.0*q0*q1*log(rxz)*uy

def repulsive_u(ep,x0,y0,z0,x1,y1,z1):
	r2 = (x1-x0)**2 + dist(y1,y0)**2 + (z1-z0)**2
	if r2 > 10.0*ep**0.1666666:
		return 0.0
	else:
		return ep/r2**6

def spring_u(h,Lmax,x0,y0,z0,x1,y1,z1):
	L = sqrt((x1-x0)**2 + dist(y1,y0)**2 + (z1-z0)**2)
	if L >= Lmax:
		return 1e10
	else:
		return -(0.5*h*Lmax**2)*log(1-(L/Lmax)**2)

def on_chain(i,N,Nl):
	return False if i < N else True

def is_endpoint(i,N,Nl):
	return True if i == N or i == N+Nl-1 or i == N+Nl or i == N+2*Nl-1 else False

def linked(i,j,N,Nl):
	if i > j:
		i, j = j, i
	if not(on_chain(i,N,Nl)) or not(on_chain(j,N,Nl)): 
		return False
	elif (i == N and j == N+Nl-1) or (i == N+Nl and j == N+2*Nl-1):
		return True
	elif abs(i-j) != 1:
		return False
	else:
		return True

def total_u(i,q0,x0,y0,z0,j,q1,x1,y1,z1,ep,h,Lmax,N,Nl):
	total = 0.0
	total += lekner_u(q0,x0,y0,z0,q1,x1,y1,z1) + repulsive_u(ep,x0,y0,z0,x1,y1,z1)
	if linked(i,j,N,Nl):
		total += spring_u(h,Lmax,x0,y0,z0,x1,y1,z1)
	return total

def delta_u(xn,x,n,ep,h,Lmax,N,Nl):
	delta_e = 0.0
	for i in range(rank,N+2*Nl,size):
		if i != n:
			new_e = total_u(n,xn[n][0], xn[n][1], xn[n][2], xn[n][3],\
							i,xn[i][0], xn[i][1], xn[i][2], xn[i][3], ep, h, Lmax, N, Nl)
			old_e = total_u(n,x[n][0], x[n][1], x[n][2], x[n][3],\
							i,x[i][0], x[i][1], x[i][2], x[i][3], ep, h, Lmax, N, Nl)
			delta_e += new_e - old_e

	return delta_e
	
def main(argv=None):

	start = time.clock()
	#Get command-line arguments
	if argv is None:
		argv = sys.argv
	
	N 	  	= int(argv[1])
	Nl		= int(argv[2])
	qci 	= float(argv[3])
	ep		= float(argv[4])
	h		= float(argv[5])
	Lmax	= float(argv[6])
	T		= int(argv[7])
	kbt0  	= float(argv[8])
	kbtf	= float(argv[9])
	R		= float(argv[10])
	fname	= argv[11]
	pOut	= int(argv[12])

	#charge on chain particles
	ql 		= -0.5*N*qci/Nl
	
	#Per-processor random seed
	random.seed(154029178285)
	random.jumpahead(rank)
	
	#Open output files
	#if gOut != 0:
	#	fg = open(fname + '.' + rank + '.growth','w')
	#if fOut != 0:
	#	ff = open(fname + '.' + rank + '.force','w')
	if pOut != 0:
		fp = open(fname + '.' + str(rank) + '.xyz','w')
	
	#Initialize particle positions
	x = zeros((N+2*Nl,4))
	for i,val in enumerate(x):
		if i < N:
			val[0] = qci
			val[1] = random.uniform(-maxR,maxR)
			val[2] = random.uniform(-Ly/2.0,Ly/2.0)
			val[3] = random.uniform(-maxR,maxR)
		elif i < N + Nl:
			val[0] = ql
			val[1] = -0.5*R
			val[2] = (i-N)*Ly/Nl - 0.5*Ly
			val[3] = 0
		else:
			val[0] = ql
			val[1] = 0.5*R
			val[2] = (i-N-Nl)*Ly/Nl - 0.5*Ly
			val[3] = 0
	
	xn = deepcopy(x)
	
	#MC loop
	acceptance_rate = 0.0;
	for s in range(T/size):
		kbt 		= kbt0#(kbtf-kbt0)*s/(T/size-1)+kbt0
		ranN 		= random.randint(0,N+2*Nl-1)
		dx, dy, dz 	= random.uniform(-1.0,1.0), random.uniform(-1.0,1.0), random.uniform(-1.0,1.0)
		step		= 0.6 if not on_chain(ranN,N,Nl) else 0.2

		xn[ranN][1] += step*dx
		xn[ranN][2] += step*dy if not is_endpoint(ranN,N,Nl) else 0
		xn[ranN][3] += step*dz
		
		DE 		= array(delta_u(xn,x,ranN,ep,h,Lmax,N,Nl), dtype='d')
		comm.Reduce([DE, MPI.DOUBLE], None, op=MPI.SUM, root=0)

		if DE < 0 or exp(-DE/kbt) > random.random():
			acceptance_rate += 1.0/T
			x[ranN][1]  = xn[ranN][1]
			x[ranN][2]  = xn[ranN][2]
			x[ranN][3]  = xn[ranN][3]
		else:
			xn[ranN][1] = x[ranN][1]
			xn[ranN][2] = x[ranN][2]
			xn[ranN][3] = x[ranN][3]
		
		if pOut != 0 and s % pOut == 0 and rank == 0:
			fp.write("%d\n" % (N + 2*Nl))
			fp.write("kbt:%f, qci:%f, ep:%f, h:%f, Lmax:%f\n" % (kbt,qci,ep,h,Lmax))
			for i,xi in enumerate(x):
				fp.write("%d %f %f %f\n" % (8 if on_chain(i,N,Nl) else 1, xi[1],xi[2],xi[3]))
	
	print acceptance_rate

	if fp != None and rank == 0:
		fp.close()
	#Close MPI

	end = time.clock()
	print '%d particles, %d steps: %.6f seconds' % (N+2*Nl,T,(end-start))
	
if __name__ == "__main__":
	main()

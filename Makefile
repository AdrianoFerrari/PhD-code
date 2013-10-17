main: 
	icpc main.c -o main.exe -limf -lgsl -lgslcblas -openmp -O3 -xHost -L/opt/intel/lib/intel64
main-debug: 
	icpc -DDEBUG -g main.c -o main.exe -limf -lgsl -lgslcblas -openmp -O3 -xHost -I${SCINET_GSL_INC} -L${SCINET_GSL_LIB}
tests:
	gcc tests.c -o tests.exe -limf -lgsl -lgslcblas -fopenmp -O3 -march=native -std=c99 -I${SCINET_GSL_INC} -L${SCINET_GSL_LIB}

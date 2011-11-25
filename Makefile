main: 
	gcc main.c -o main.exe -limf -lgsl -lgslcblas -fopenmp -O4 -std=c99 -I${SCINET_GSL_INC} -L${SCINET_GSL_LIB}
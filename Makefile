main: 
	gcc -g main.c -o main -lm -lgsl -lgslcblas -fopenmp -O4 -std=c99 -I${SCINET_GSL_INC} -L${SCINET_GSL_LIB}
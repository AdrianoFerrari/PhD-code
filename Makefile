main: 
	gcc main.c -o main.exe -limf -lgsl -lgslcblas -fopenmp -O4 -std=c99 -I${SCINET_GSL_INC} -L${SCINET_GSL_LIB}
force_calc:
	gcc force_calc.c -o force_calc.exe -limf -lgsl -lgslcblas -fopenmp -O4 -std=c99 -I${SCINET_GSL_INC} -L${SCINET_GSL_LIB}
tests:
	gcc tests.c -o tests.exe -limf -lgsl -lgslcblas -fopenmp -O4 -std=c99 -I${SCINET_GSL_INC} -L${SCINET_GSL_LIB}

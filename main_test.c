#include "main.c"

#include <assert.h>
#include <stdbool.h>

static void test_dummy_function() {
  assert(dummy_function() == 1.234);
}

static void test_lekner_potential_1() {
  assert(lekner_potential( 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0 ) >= -0.0142195*1.005
	 &&
	 lekner_potential( 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0 ) <= -0.0142195*0.995
	 );
}
static void test_lekner_potential_2() {
  assert(lekner_potential( 1.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0, 1.0 ) >= 0.105726*0.995
	 &&
	 lekner_potential( 1.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0, 1.0 ) <= 0.105726*1.005
	 );
}
static void test_repulsive_potential_1() {
  assert(repulsive_potential(1.23, 4.21, -5.0, 3.22, 1.044, 0.0, 23) == 0.0);
}
static void test_repulsive_potential_2() {
  assert(repulsive_potential(1.23,4.21,-1.0,3.22,1.044,0.0,23) >= 0.00000203246*0.995
	 &&
	 repulsive_potential(1.23,4.21,-1.0,3.22,1.044,0.0,23) <= 0.00000203246*1.005
	 );
}
static void test_repulsive_potential_3() {
  assert(repulsive_potential(0.35,0.90,0.70,-0.27,-0.07,0.49,0.81) >= 0.12283*0.995
	 &&
	 repulsive_potential(0.35,0.90,0.70,-0.27,-0.07,0.49,0.81) <= 0.12283*1.005
	 );
}
static void test_repulsive_potential_4() {
  assert(repulsive_potential(1.0,1.0,1.0,1.0,1.0,1.0,1.0) == 10e20);
}
static void test_kiss_rng() {
  int i; 
  for(i=0;i<100000000;i++) t=KISS;
  assert(t==1666297717051644203ULL);
}

int main() {
  //test_dummy_function();
  //test_lekner_potential_1();
  //test_lekner_potential_2();
  //test_repulsive_potential_1();
  //test_repulsive_potential_2();  
  //test_repulsive_potential_3();
  //test_repulsive_potential_4();
  //test_kiss_rng();
}

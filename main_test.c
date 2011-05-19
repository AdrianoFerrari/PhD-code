#include "main.c"

#include <assert.h>
#include <stdbool.h>

static void test_dummy_function() {
  assert(dummy_function() == 1.234);
}

static void test_lekner_potential_1() {
  double val = lekner_potential( 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0 );
  double expected = -0.0142195;
  assert(val >= expected*1.005 && val <= expected*0.995);
}
static void test_lekner_potential_2() {
  double val = lekner_potential( 1.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0, 1.0 );
  double expected = 0.105726;
  assert(val >= expected*0.995 && val <= expected*1.005);
}
static void test_repulsive_potential_1() {
  assert(repulsive_potential(1.23, 4.21, -5.0, 3.22, 1.044, 0.0, 23) == 0.0);
}
static void test_repulsive_potential_2() {
  double val = repulsive_potential(1.23,4.21,-1.0,3.22,1.044,0.0,23);
  double expected = 0.00000203246;
  assert(val >= expected*0.995 && val <= expected*1.005);
}
static void test_repulsive_potential_3() {
  double val = repulsive_potential(0.35,0.90,0.70,-0.27,-0.07,0.49,0.81);
  double expected = 0.12283;
  assert(val >= expected*0.995 && val <= expected*1.005);
}
static void test_repulsive_potential_4() {
  assert(repulsive_potential(1.0,1.0,1.0,1.0,1.0,1.0,1.0) == 10e20);
}
static void test_kiss_rng() {
  int i; 
  for(i=0;i<100000000;i++) t=KISS;
  assert(t==1666297717051644203ULL);
}
static void test_kiss_rng_seed() {
  int i; 
  for(i=0;i<1000;i++) t=KISS;  
  settable(1234567890987654321ULL, 123456123456123456ULL, 362436362436362436ULL, 1066149217761810ULL);
  for(i=0;i<100000000;i++) t=KISS;
  assert(t==1666297717051644203ULL);
}
static void test_pair_potential_energy_1() {
  double val = pair_potential_energy(-0.043, 0.256, 0.733, -0.209, 0.957, -0.555, 0.916, 0.637, 0.208);
  double expected = 0.00474872;
  assert(val >= expected*0.995
	 &&
	 val <= expected*1.005
	 );
}
static void test_pair_potential_energy_2() {
  double val = pair_potential_energy(-0.668, 0.137, 0.260, 0.511, -0.531, 0.082, -0.512, 0.281, 0.410);
  double expected = 5.90082;
  assert(val >= expected*0.995
	 &&
	 val <= expected*1.005
	 );
}
static void test_single_particle_energy_1() {
  int size = 6;
  double q[6] = {-0.485249, 0.492407, 0.667379, -0.845798, -0.218984, 0.298472};
  double x[6][3] = {{22.7278, 20.5178, -4.34251}, {2.05805, 6.01785, 12.1014}, {23.2642, 21.6427, -31.9224}, {27.7179, 21.0241, 10.711}, {18.0486, 9.40886, 8.45633}, {-30.2842, 11.2869, 11.9472}};
  double eps = 1.234;
  double val = particle_total_potential(x,q,size,eps,4);
  double expected = 0.0308394;
  assert(val >= expected*0.995
	 &&
	 val <= expected*1.005
	 );
}
static void test_single_particle_energy_2() {
  int size = 6;
  double q[6] = {0.161532, 0.395393, 0.59779, -0.0170533, -0.536349, -0.7383};
  double x[6][3] = {{15.8407, 17.7708, -20.4505}, {19.3111, 6.65337, -28.5583}, {19.5652, 12.5071, 9.00004}, {-3.57101, 19.715, 5.84241}, {22.9855,  16.549, -14.7609}, {-19.6095, 7.88197, 7.93345}};
  double eps = 2.107;
  double val = particle_total_potential(x,q,size,eps,1);
  double expected = 0.062502;
  assert(val >= expected*0.995
	 &&
	 val <= expected*1.005
	 );
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
  test_kiss_rng_seed();
  //test_pair_potential_energy_1();  
  //test_pair_potential_energy_2();
  //test_single_particle_energy_1();  
  //test_single_particle_energy_2();
}

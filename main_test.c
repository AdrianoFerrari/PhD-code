#include "main.c"

#include <assert.h>

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
static void test_line_potential_1() {
  double val = line_potential(1.0,2.0,3.0,4.0, 2.0,0.0,0.0, 0.5, 2.0);
  double expected = -5.99146;
  assert(val <= expected*0.995
	 &&
	 val >= expected*1.005
	 );
}
static void test_line_potential_2() {
  double val = line_potential(-0.223613, -0.613204, -0.500095, -0.811417, -0.712559, -0.0448241, -0.984596, 0.856894, 0.609515);
  double expected = 0.164019;
  assert(val >= expected*0.995
	 &&
	 val <= expected*1.005
	 );
}
static void test_line_repulsive_potential_1() {
  double val = line_repulsive_potential(0.664857, 0.302964, -0.762206, -0.842198, 0.262977, -0.21951, 0.843407, 0.353052);
  double expected = 0.911892; 
  assert(val >= expected*0.995
	 &&
	 val <= expected*1.005
	 );
}
static void test_line_repulsive_potential_2() {
  double val = line_repulsive_potential(0.128893,1.3044,-11.3497,-4.12898,-11.591,1.11956,10.933,8.8148);
  double expected = 0.0; 
  assert(val >= expected*0.995
	 &&
	 val <= expected*1.005
	 );
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
static void test_kiss_uni() {
  double sum = 0.0;
  int i;
  for(i=0; i<100000;i++){
    sum += UNI;
  }
  assert(sum/100000 >= 0.5*0.995
	 &&
	 sum/100000 <= 0.5*1.005
	 );
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
static void test_particle_total_pair_potential_1() {
  int size = 6;
  double q[6] = {-0.485249, 0.492407, 0.667379, -0.845798, -0.218984, 0.298472};
  double x[6][3] = {{22.7278, 20.5178, -4.34251}, {2.05805, 6.01785, 12.1014}, {23.2642, 21.6427, -31.9224}, {27.7179, 21.0241, 10.711}, {18.0486, 9.40886, 8.45633}, {-30.2842, 11.2869, 11.9472}};
  double eps = 1.234;
  double qL = 1.00; double xL = 1.234; double zL = -0.324; double A = 0.12; double lambda = 3.5;
  double val = particle_total_pair_potential(x,q,size,eps,1);
  double expected = 0.053867;
  assert(val >= expected*0.995
	 &&
	 val <= expected*1.005
	 );
}
static void test_particle_total_potential_1() {
  int size = 6;
  double q[6] = {-0.485249, 0.492407, 0.667379, -0.845798, -0.218984, 0.298472};
  double x[6][3] = {{22.7278, 20.5178, -4.34251}, {2.05805, 6.01785, 12.1014}, {23.2642, 21.6427, -31.9224}, {27.7179, 21.0241, 10.711}, {18.0486, 9.40886, 8.45633}, {-30.2842, 11.2869, 11.9472}};
  double eps = 1.234;
  double qL = 1.00; double xL1 = 1.234; double zL1 = -0.324; double A1 = 0.12; double lambda = 3.5;
  double xL2 = -1.234; double zL2 = 0.324; double A2 = 0.12;
  double val = particle_total_potential(x,q,size,1,eps,qL,xL1,zL1,xL2,zL2,A1,A2,lambda);
  double expected = -4.89555;
  assert(val <= expected*0.995
	 &&
	 val >= expected*1.005
	 );
}
static void test_particle_total_potential_2() {
  int size = 6;
  double q[6] = {0.161532, 0.395393, 0.59779, -0.0170533, -0.536349, -0.7383};
  double x[6][3] = {{15.8407, 17.7708, -20.4505}, {19.3111, 6.65337, -28.5583}, {19.5652, 12.5071, 9.00004}, {-3.57101, 19.715, 5.84241}, {22.9855,  16.549, -14.7609}, {-19.6095, 7.88197, 7.93345}};
  double eps = 2.107;
  double qL = 1.20; double xL1 = 1.234; double zL1 = -0.324; double A1 = 0.12; double lambda = 3.5;
  double xL2 = -1.234; double zL2 = 0.324; double A2 = 0.12;
  double val = particle_total_potential(x,q,size,1,eps,qL,xL1,zL1,A1,xL2,zL2,A2,lambda);
  double expected = -6.65627;
  assert(val <= expected*0.995
	 &&
	 val >= expected*1.005
	 );
}
static void test_energy_difference_1() {
  int size = 6;
  double q[6] = {0.161532, 0.395393, 0.59779, -0.0170533, -0.536349, -0.7383};
  double x[6][3] = {{15.8407, 17.7708, -20.4505}, {19.3111, 6.65337, -28.5583}, {19.5652, 12.5071, 9.00004}, {-3.57101, 19.715, 5.84241}, {22.9855,  16.549, -14.7609}, {-19.6095, 7.88197, 7.93345}};
  double xn[6][3] = {{15.8407, 17.7708, -20.4505}, {19.3111, 6.65337, -28.5583}, {19.5652, 12.5071, 9.00004}, {31.0154, 5.28353, 14.9521}, {22.9855,  16.549, -14.7609}, {-19.6095, 7.88197, 7.93345}};
  double eps = 2.107;
  double qL = 1.20; double xL1 = 1.234; double zL1 = -0.324; double A1 = 0.12; double lambda = 3.5;
  double xL2 = -1.234; double zL2 = 0.324; double A2 = 0.12;
  double val = energy_difference(x,xn,q,size,3,eps,qL,xL1,zL1,A1,xL2,zL2,A2,lambda);
  double expected = 0.130702;
  assert(val >= expected*0.995
	 &&
	 val <= expected*1.005
	 );
}
static void test_energy_difference_2() {
  int size = 6;
  double q[6] = {0.161532, 0.395393, 0.59779, -0.0170533, -0.536349, -0.7383};
  double x[6][3] = {{15.8407, 17.7708, -20.4505}, {19.3111, 6.65337, -28.5583}, {19.5652, 12.5071, 9.00004}, {-3.57101, 19.715, 5.84241}, {22.9855,  16.549, -14.7609}, {-19.6095, 7.88197, 7.93345}};
  double xn[6][3] = {{15.8407, 17.7708, -20.4505}, {19.3111, 6.65337, -28.5583}, {19.5652, 12.5071, 9.00004}, {19.3, 6.6, -28.0}, {22.9855,  16.549, -14.7609}, {-19.6095, 7.88197, 7.93345}};
  double eps = 2.107;
  double qL = 1.20; double xL1 = 1.234; double zL1 = -0.324; double A1 = 0.12; double lambda = 3.5;
  double xL2 = -1.234; double zL2 = 0.324; double A2 = 0.12;
  double val = energy_difference(x,xn,q,size,3,eps,qL,xL1,zL1,A1,xL2,zL2,A2,lambda);
  double expected = 2170.47;
  assert(val >= expected*0.995
	 &&
	 val <= expected*1.005
	 );
}
static void test_go_ahead_1() {
  double kbT = 1.0;
  double De = -0.123;
  bool val = go_ahead(De,kbT);
  assert( val == true);
}
static void test_go_ahead_2() {
  double kbT = 1.0;
  double De = 1e9;
  bool val = go_ahead(De,kbT);
  assert( val == false);
}
static void test_go_ahead_3() {
  double sum = 0.0;
  double expected = 0.2236;
  int i;
  for(i=0;i<100000;i++) {
    sum += go_ahead(1.5,1.0);
  }
  assert(sum/100000 >= expected*0.95
	 &&
	 sum/100000 <= expected*1.05
	 );
}
static void test_go_ahead_4() {
  double sum = 0.0;
  double expected = 0.011193;
  int i;
  for(i=0;i<100000;i++) {
    sum += go_ahead(4.5,1.0);
  }
  assert(sum/100000 >= expected*0.95
	 &&
	 sum/100000 <= expected*1.05
	 );
}
static void test_line_growth_at_pt_from_particle_1() {
  double s = 1.2; int Ns = 14;
  double q = 1.0;
  double x = 1.234; double y = 1.234; double z = -1.234;
  double qL = -1.2; double xL = -0.5; double zL = 0.0; double A = 1.2; double lambda = 3.5;
  double val = line_growth_at_pt_from_particle(1.2, 1, 1.234, 1.234, -1.234, -1.2, -0.5, 0.0, 1.2, 3.5, 14);
  double expected = -0.014285;
  assert(val <= expected*0.995
	 &&
	 val >= expected*1.005
	 );
}
static void  test_xforce_on_seg_due_to_part_1() {
  double s = 0.22;
  double q = 1.0; double x = 2.1; double y = 1.5; double z = 0.0;
  double qL = -16.0; double xL = -1.0; double zL = 0.0;
  double ep = 1.2; double A = 0.1; double lambda = 3.5; int Ns = 23; int M = 7;
  double val = xforce_on_seg_due_to_part(s,q,x,y,z,qL,xL,zL,ep,A,lambda,Ns,M);
  double expected = 0.455629;
  assert(val >= expected*0.995
	 &&
	 val <= expected*1.005
	 );
}
static void  test_xforce_on_seg_due_to_part_2() {
  double s = 0.22;
  double q = 1.0; double x = -0.9; double y = 1.5; double z = 0.0;
  double qL = -16.0; double xL = -1.0; double zL = 0.0;
  double ep = 1.2; double A = 0.1; double lambda = 3.5; int Ns = 23; int M = 7;
  double val = xforce_on_seg_due_to_part(s,q,x,y,z,qL,xL,zL,ep,A,lambda,Ns,M);
  double expected = -7.2836e13;
  assert(val <= expected*0.995
	 &&
	 val >= expected*1.005
	 );
}
static void  test_xforce_on_seg_due_to_part_3() {
  double s = 0.22;
  double q = 1.0; double x = -1.1; double y = 1.5; double z = 0.0;
  double qL = -16.0; double xL = -1.0; double zL = 0.0;
  double ep = 1.2; double A = 0.1; double lambda = 3.5; int Ns = 23; int M = 7;
  double val = xforce_on_seg_due_to_part(s,q,x,y,z,qL,xL,zL,ep,A,lambda,Ns,M);
  double expected = 2.95639e14;
  assert(val >= expected*0.995
	 &&
	 val <= expected*1.005
	 );
}
static void  test_xforce_on_seg_due_to_part_4() {
  double s = 0.22;
  double q = 1.0; double x = -2.9; double y = 1.5; double z = 0.0;
  double qL = -16.0; double xL = -1.0; double zL = 0.0;
  double ep = 1.2; double A = 0.1; double lambda = 3.5; int Ns = 23; int M = 7;
  double val = xforce_on_seg_due_to_part(s,q,x,y,z,qL,xL,zL,ep,A,lambda,Ns,M);
  double expected = -0.482637;
  assert(val <= expected*0.995
	 &&
	 val >= expected*1.005
	 );
}
static void  test_xforce_on_seg_due_to_line_1() {
  double s = 0.63;
  double qL = -16.0; double xL1 = 1.5; double zL1 = 0.0;
  double xL2 = -1.5; double zL2 = 0.0;
  double A = 0.1; double lambda = 3.5; int Ns = 23;
  double val = xforce_on_seg_due_to_line(s,qL,xL1,zL1,xL2,zL2,A,lambda,Ns);
  double expected = 172.873;
  assert(val >= expected*0.995
	 &&
	 val <= expected*1.005
	 );
}
static void  test_xforce_on_seg_due_to_line_2() {
  double s = 0.63;
  double qL = -16.0; double xL1 = 0.9; double zL1 = 0.0;
  double xL2 = -0.9; double zL2 = 0.0;
  double A = 0.1; double lambda = 3.5; int Ns = 23;
  double val = xforce_on_seg_due_to_line(s,qL,xL1,zL1,xL2,zL2,A,lambda,Ns);
  double expected = 282.606;
  assert(val >= expected*0.995
	 &&
	 val <= expected*1.005
	 );
}
static void  test_xforce_on_seg_due_to_line_3() {
  double s = 0.63;
  double qL = -16.0; double xL1 = -3.9; double zL1 = 0.0;
  double xL2 = 3.9; double zL2 = 0.0;
  double A = 0.1; double lambda = 3.5; int Ns = 23;
  double val = xforce_on_seg_due_to_line(s,qL,xL1,zL1,xL2,zL2,A,lambda,Ns);
  double expected = -69.2989;
  assert(val <= expected*0.995
	 &&
	 val >= expected*1.005
	 );
}
static void test_growth_on_seg_due_to_part_1() {
  double s = 0.22;
  double q = 1.0; double x = 2.1; double y = 1.5; double z = 0.0;
  double qL = -16.0; double xL = -1.0; double zL = 0.0;
  double ep = 1.2; double A = 0.1; double lambda = 3.5; int Ns = 23; int M = 7;
  double val = growth_on_seg_due_to_part(s,q,x,y,z,qL,xL,zL,ep,A,lambda,Ns,M);
  double expected = 0.0336975;
  assert(val >= expected*0.995
	 &&
	 val <= expected*1.005
	 );  
}
static void test_growth_on_seg_due_to_part_2() {
  double s = 0.22;
  double q = 1.0; double x = -0.9; double y = 1.5; double z = 0.0;
  double qL = -16.0; double xL = -1.0; double zL = 0.0;
  double ep = 1.2; double A = 0.1; double lambda = 3.5; int Ns = 23; int M = 7;
  double val = growth_on_seg_due_to_part(s,q,x,y,z,qL,xL,zL,ep,A,lambda,Ns,M);
  double expected = -1.32202e16;
  assert(val <= expected*0.995
	 &&
	 val >= expected*1.005
	 );  
}
static void test_growth_on_seg_due_to_part_3() {
  double s = 0.22;
  double q = 1.0; double x = -1.1; double y = 1.5; double z = 0.0;
  double qL = -16.0; double xL = -1.0; double zL = 0.0;
  double ep = 1.2; double A = 0.1; double lambda = 3.5; int Ns = 23; int M = 7;
  double val = growth_on_seg_due_to_part(s,q,x,y,z,qL,xL,zL,ep,A,lambda,Ns,M);
  double expected = -2.81701e16;
  assert(val <= expected*0.995
	 &&
	 val >= expected*1.005
	 );  
}
static void test_growth_on_seg_due_to_part_4() {
  double s = 0.22;
  double q = 1.0; double x = -2.9; double y = 1.5; double z = 0.0;
  double qL = -16.0; double xL = -1.0; double zL = 0.0;
  double ep = 1.2; double A = 0.1; double lambda = 3.5; int Ns = 23; int M = 7;
  double val = growth_on_seg_due_to_part(s,q,x,y,z,qL,xL,zL,ep,A,lambda,Ns,M);
  double expected = 0.0435877;
  assert(val >= expected*0.995
	 &&
	 val <= expected*1.005
	 );  
}
static void test_growth_on_seg_due_to_line_1() {
  double s = 0.63;
  double qL = -16.0; double xL1 = 1.5; double zL1 = 0.0;
  double xL2 = -1.5; double zL2 = 0.0;
  double A = 0.1; double lambda = 3.5; int Ns = 23;
  double val = growth_on_seg_due_to_line(s,qL,xL1,zL1,xL2,zL2,A,lambda,Ns);
  double expected = -57.6243;
  assert(val <= expected*0.995
	 &&
	 val >= expected*1.005
	 );
}
static void test_growth_on_seg_due_to_line_2() {
  double s = 0.63;
  double qL = -16.0; double xL1 = 0.9; double zL1 = 0.0;
  double xL2 = -0.9; double zL2 = 0.0;
  double A = 0.1; double lambda = 3.5; int Ns = 23;
  double val = growth_on_seg_due_to_line(s,qL,xL1,zL1,xL2,zL2,A,lambda,Ns);
  double expected = -157.003;
  assert(val <= expected*0.995
	 &&
	 val >= expected*1.005
	 );
}
static void test_growth_on_seg_due_to_line_3() {
  double s = 0.63;
  double qL = -16.0; double xL1 = -3.9; double zL1 = 0.0;
  double xL2 = 3.9; double zL2 = 0.0;
  double A = 0.1; double lambda = 3.5; int Ns = 23;
  double val = growth_on_seg_due_to_line(s,qL,xL1,zL1,xL2,zL2,A,lambda,Ns);
  double expected = -8.88447;
  assert(val <= expected*0.995
	 &&
	 val >= expected*1.005
	 );
}
static void test_growth_on_seg_1() {
  double s = 0.43;
  double q[6] = {0.161532, 0.395393, 0.59779, 0.0170533, 0.536349, 0.7383};
  double x[6][3] = {{15.8407, 17.7708, -20.4505}, {19.3111, 6.65337, -28.5583}, {19.5652,12.5071, 9.00004}, {-3.57101, 19.715, 5.84241}, {22.9855, 16.549, -14.7609}, {-19.6095, 7.88197, 7.93345}};
  double ep = 2.107; double qL = -16.0; double xL1 = 1.234; double zL1 = 0.0; double A = 0.1; double lambda = 3.5;
  double xL2 = -1.234; double zL2 = 0.0; int Ns = 23; int M = 7;
  double val = growth_on_seg(s,q,x,6,qL,xL1,zL1,xL2,zL2,ep,A,lambda,Ns,M);
  double expected = -88.853;
  assert(val <= expected*0.995
	 &&
	 val >= expected*1.005
	 );
}
static void test_growth_on_seg_2() {
  double s = 0.43;
  double q[6] = {0.161532, 0.395393, 0.59779, 0.0170533, 0.536349, 0.7383};
  double x[6][3] = {{1.717, -0.619, 1.483}, {0.28, -0.341, -2.115}, {0.908, -4.512, 2.473}, {-0.913, 3.142, 3.084}, {1.731, -4.153, -0.95}, {2.642, 3.118, -1.636}};
  double ep = 2.107; double qL = -16.0; double xL1 = 25; double zL1 = 0.0; double A = 0.1; double lambda = 3.5;
  double xL2 = -25; double zL2 = 0.0; int Ns = 23; int M = 7;
  double val = growth_on_seg(s,q,x,6,qL,xL1,zL1,xL2,zL2,ep,A,lambda,Ns,M);
  double expected = -0.207639;
  assert(val <= expected*0.995
	 &&
	 val >= expected*1.005
	 );
}
static void test_growth_rate_1() {
  double q[6] = {0.161532, 0.395393, 0.59779, 0.0170533, 0.536349, 0.7383};
  double x[6][3] = {{15.8407, 17.7708, -20.4505}, {19.3111, 6.65337, -28.5583}, {19.5652,12.5071, 9.00004}, {-3.57101, 19.715, 5.84241}, {22.9855,16.549, -14.7609}, {-19.6095, 7.88197, 7.93345}};
  double ep = 2.107; double qL = -16.0; double xL1 = 1.234; double zL1 = 0.0; double A = 0.1; double lambda = 3.5;
  double xL2 = -1.234; double zL2 = 0.0; int Ns = 23; int M = 7;
  double val = growth_rate(q,x,6,qL,xL1,zL1,xL2,zL2,ep,A,lambda,Ns,M);
  double expected = -87.5393;
  assert(val <= expected*0.995
	 &&
	 val >= expected*1.005
	 );
}
static void test_growth_rate_2() {
  double q[6] = {0.161532, 0.395393, 0.59779, 0.0170533, 0.536349, 0.7383};
  double x[6][3] = {{1.717, -0.619, 1.483}, {0.28, -0.341, -2.115}, {0.908, -4.512, 2.473}, {-0.913, 3.142, 3.084}, {1.731, -4.153, -0.95}, {2.642, 3.118, -1.636}};
  double ep = 2.107; double qL = -16.0; double xL1 = 25; double zL1 = 0.0; double A = 0.1; double lambda = 3.5;
  double xL2 = -25; double zL2 = 0.0; int Ns = 23; int M = 7;
  double val = growth_rate(q,x,6,qL,xL1,zL1,xL2,zL2,ep,A,lambda,Ns,M);
  double expected = -0.207572;
  assert(val <= expected*0.995
	 &&
	 val >= expected*1.005
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
  //test_line_potential_1();
  //test_line_potential_2();
  //test_line_repulsive_potential_1();
  //test_line_repulsive_potential_2();
  //test_kiss_rng();
  //test_kiss_rng_seed();
  //test_kiss_uni();
  //test_pair_potential_energy_1();  
  //test_pair_potential_energy_2();
  //test_particle_total_pair_potential_1();
  //test_particle_total_potential_1();  
  //test_particle_total_potential_2();
  //test_energy_difference_1();
  //test_energy_difference_2();
  //test_go_ahead_1();
  //test_go_ahead_2();
  //test_go_ahead_3();
  //test_go_ahead_4();
  //test_xforce_on_seg_due_to_part_1();
  //test_xforce_on_seg_due_to_part_2();
  //test_xforce_on_seg_due_to_part_3();
  //test_xforce_on_seg_due_to_part_4();
  //test_xforce_on_seg_due_to_line_1();
  //test_xforce_on_seg_due_to_line_2();
  //test_xforce_on_seg_due_to_line_3();
  //test_growth_on_seg_due_to_part_1();
  //test_growth_on_seg_due_to_part_2();
  //test_growth_on_seg_due_to_part_3();
  //test_growth_on_seg_due_to_part_4();
  //test_growth_on_seg_due_to_line_1();
  test_growth_on_seg_due_to_line_2();
  test_growth_on_seg_due_to_line_3();
  test_growth_on_seg_1();
  test_growth_on_seg_2();
  test_growth_rate_1();
  test_growth_rate_2();
}

#include "functions.c"
#include <assert.h>

static void test_linked() {
  assert(linked(6, 7, 4, 5));
  assert(linked(7, 6, 4, 5));
  assert(linked(4, 8, 4, 5));
  assert(linked(9,13, 4, 5));
  assert(linked(9,10, 4, 5));

  assert(!linked(6, 8, 4, 5));
  assert(!linked(6, 4, 4, 5));
  assert(!linked(9,11, 4, 5));
}
static void test_dist() {
  assert(dist(0.0, 1.0) == 1.0);
  assert(dist(0.0, 3.0) == 3.0);
  assert(dist(5.0, 6.0) == 1.0);
  assert(dist(-5.0, 5.0) == 10.0);
  assert(dist(-5.0, 6.0) == 11.0);
  assert(dist(-5.0, 7.0) == 12.0);
  assert(dist(-5.0, 8.0) == 11.0);
  assert(dist(0.0, 13.0) == 11.0);
}
static void test_on_chain() {
  assert(!on_chain(0,32));
  assert(!on_chain(0,32));

  assert( on_chain(32,32));
  assert( on_chain(35,32));
} 
static void test_is_endpoint() {
  assert(!is_endpoint(0, 6, 5));
  assert(!is_endpoint(5, 6, 5));
  assert(!is_endpoint(7, 6, 5));
  assert(!is_endpoint(12, 6, 5));
  assert(!is_endpoint(14, 6, 5));

  assert( is_endpoint(6, 6, 5));
  assert( is_endpoint(10, 6, 5));
  assert( is_endpoint(11, 6, 5));
  assert( is_endpoint(15, 6, 5));
} 
static void test_theta_du() {
  double val = theta_du(1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0);
  double exp = 2.4674;
  double tol = 0.001;
  assert(val >= exp-tol && val <= exp+tol);

  val = theta_du(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 23.0, 0.0, 0.0, 1.0, 0.0);
  exp = 2.4674;
  tol = 0.001;
  assert(val >= exp-tol && val <= exp+tol);
} 
static void test_delta_u() {
  int N = 2; int Nl = 2;
  double *q; 
  double **x; double **xn;
  q = malloc((N+2*Nl)*sizeof(double));
  x = malloc((N+2*Nl)*sizeof(double *));
  xn = malloc((N+2*Nl)*sizeof(double *));

  for(int i=0;i<N+2*Nl;i++) {
    x[i] = malloc(3*sizeof(double));
    xn[i] = malloc(3*sizeof(double));
  }

  q[0] = q[1] = 1.0;
  q[2] = q[3] = q[4] = q[5] = -0.5;

  x[0][0] = xn[0][0] = 0.0; x[0][1] = xn[0][1] = 1.0; x[0][2] = xn[0][2] = 2.0;
  x[1][0] = xn[1][0] = 1.0; x[1][1] = xn[1][1] = 2.0; x[1][2] = xn[1][2] = 3.0;
  x[2][0] = xn[2][0] = 2.0; x[2][1] = xn[2][1] = 3.0; x[2][2] = xn[2][2] = 4.0;
  x[3][0] = xn[3][0] = 3.0; x[3][1] = xn[3][1] = 4.0; x[3][2] = xn[3][2] = 5.0;
  x[4][0] = xn[4][0] = 4.0; x[4][1] = xn[4][1] = 5.0; x[4][2] = xn[4][2] = 6.0;
  x[5][0] = xn[5][0] = 5.0; x[5][1] = xn[5][1] = 6.0; x[5][2] = xn[5][2] = 7.0;

  xn[0][2] += 0.1;

  double val = delta_u(x,xn,q,0,1.0,1.0,1.0,1.0,2,2);
  double exp = 0.0180702;
  double tol = 0.00001;
  assert(val >= exp-tol && val <= exp+tol);

  x[0][0] = xn[0][0] = 0.22301800; x[0][1] = xn[0][1] = 0.954541000; x[0][2] = xn[0][2] = 0.726033;
  x[1][0] = xn[1][0] = 0.00665308; x[1][1] = xn[1][1] = 0.138399000; x[1][2] = xn[1][2] = 0.332936;
  x[2][0] = xn[2][0] = 0.62877400; x[2][1] = xn[2][1] = 0.000933274; x[2][2] = xn[2][2] = 0.344204;
  x[3][0] = xn[3][0] = 0.70361800; x[3][1] = xn[3][1] = 0.788098000; x[3][2] = xn[3][2] = 0.922483;
  x[4][0] = xn[4][0] = 0.55393200; x[4][1] = xn[4][1] = 0.201894000; x[4][2] = xn[4][2] = 0.963744;
  x[5][0] = xn[5][0] = 0.90615800; x[5][1] = xn[5][1] = 0.305324000; x[5][2] = xn[5][2] = 0.904319;

  xn[0][0] += 0.2;

  val = delta_u(x,xn,q,0,3.0,4.0,5.0,1.0,2,2);
  exp = 318057;
  tol = 7.0;
  assert(val >= exp-tol && val <= exp+tol);

  x[0][0] = xn[0][0] = 0.0; x[0][1] = xn[0][1] = 1.0; x[0][2] = xn[0][2] = 2.0;
  x[1][0] = xn[1][0] = 1.0; x[1][1] = xn[1][1] = 2.0; x[1][2] = xn[1][2] = 3.0;
  x[2][0] = xn[2][0] = 2.0; x[2][1] = xn[2][1] = 3.0; x[2][2] = xn[2][2] = 4.0;
  x[3][0] = xn[3][0] = 3.0; x[3][1] = xn[3][1] = 4.0; x[3][2] = xn[3][2] = 5.0;
  x[4][0] = xn[4][0] = 4.0; x[4][1] = xn[4][1] = 5.0; x[4][2] = xn[4][2] = 6.0;
  x[5][0] = xn[5][0] = 5.0; x[5][1] = xn[5][1] = 6.0; x[5][2] = xn[5][2] = 7.0;

  xn[1][0] -= 0.9;
  xn[1][1] -= 0.9;
  xn[1][2] -= 0.9;

  val = delta_u(x,xn,q,1,1.0,1.0,1.0,1.0,2,2);
  exp = 1.37174e9;
  tol = 0.001e9;
  assert(val >= exp-tol && val <= exp+tol);
} 
static void test_ran_xz() {
  double r;
  for(int i = 0; i < 3000; i++) {
    r = ran_xz(24.0); 
    assert(r>=-24.0 && r <=24.0);
  }
}
static void test_ran_y() {
  double r;
  for(int i = 0; i < 3000; i++) {
    r = ran_y(24.0); 
    assert(r>=0.0 && r <=24.0);
  }
}
static void test_ran_particle() {
  double r;
  for(int i = 0; i < 3000; i++) {
    r = ran_particle(24); 
    assert(r>=0 && r <=23);
  }
}
static void test_ran_u() {
  double r;
  double sum;
  for(int i = 0; i < 3000; i++) {
    r = ran_u(); 
    assert(r>=0.0 && r <=1.0);
    sum += r;
  }
  assert(sum/3000.0 >= 0.49 && sum/3000.0 <= 0.51);
}
static void test_ran_du() {
  double r;
  double sum;
  for(int i = 0; i < 3000; i++) {
    r = ran_du(); 
    assert(r>= -1.0 && r <=1.0);
    sum += r;
  }
  assert(sum/3000.0 >= -0.017 && sum/3000.0 <= 0.017);
}
static void test_lekner_fx(){
  double val = lekner_fx(1.0, -1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0); 
  double exp = -0.241112;
  double tol = 0.001;
  assert(val >= exp-tol && val <= exp+tol);

  val = lekner_fx(-146.44, 80.54, -43.32, 755.8, -839.1, 388.14, 965.54, 712.41); 
  exp = -32.6398;
  tol = 0.001;
  assert(val >= exp-tol && val <= exp+tol);

  val = lekner_fx(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0); 
  assert(isnan(val));

  val = lekner_fx(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 3.0, 0.0); 
  exp = 0.0;
  tol = 0.001;
  assert(val >= exp-tol && val <= exp+tol);
}
static void test_rep_fx() {
  double val = rep_fx(1.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0); 
  double exp = 0.00146484;
  double tol = 0.00001;
  assert(val >= exp-tol && val <= exp+tol);

  val = rep_fx(1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0); 
  exp = 0.0;
  tol = 0.001;
  assert(val >= exp-tol && val <= exp+tol);
}
static void test_spring_fx() {
  double val = spring_fx(1.0, 1.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0); 
  double exp = 0.0;
  double tol = 0.001;
  assert(val >= exp-tol && val <= exp+tol);

  val = spring_fx(1.0, 2.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0); 
  exp = 0.0;
  tol = 0.001;
  assert(val >= exp-tol && val <= exp+tol);

  val = spring_fx(1.0, 2.0, -0.9, 0.0, 0.0, 0.9, 0.0, 0.0); 
  exp = 9.47368; 
  tol = 0.001;
  assert(val >= exp-tol && val <= exp+tol);
}
int main(int argc, char **argv) {
  test_linked();
  test_dist();
  test_on_chain();
  test_is_endpoint();
  test_theta_du();
  test_delta_u();
  test_ran_xz();
  test_ran_y();
  test_ran_particle();
  test_ran_u();
  test_ran_du();

  test_lekner_fx();
  test_rep_fx();
  test_spring_fx();

  return 1;
}

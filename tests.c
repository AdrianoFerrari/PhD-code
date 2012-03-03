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
/*  double val = delta_u(1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0);
  double exp = 2.4674;
  double tol = 0.001;
  assert(val >= exp-tol && val <= exp+tol);

  val = delta_u(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 23.0, 0.0, 0.0, 1.0, 0.0);
  exp = 2.4674;
  tol = 0.001;
  assert(val >= exp-tol && val <= exp+tol);*/
  printf("Untested: delta_u\n");
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

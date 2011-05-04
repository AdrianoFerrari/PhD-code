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
int main() {
  test_dummy_function();
  test_lekner_potential_1();
  test_lekner_potential_2();
}

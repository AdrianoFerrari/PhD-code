#include "main.c"

#include <assert.h>
#include <stdbool.h>

static void test_dummy_function() {
  assert(dummy_function() == 1.234);
}

static void test_lekner_potential() {
  assert( lekner_potential( 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 ) == 1.234 );
}

int main() {
  test_dummy_function();
}

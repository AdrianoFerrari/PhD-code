#include "main.h"

#include <assert.h>
#include <stdbool.h>

static void test_lekner_potential() {
    assert( lekner_potential( 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 ) == 1.234 );
}

int main() {
    test_lekner_potential();
}

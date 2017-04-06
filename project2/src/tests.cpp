#include "tests.h" // header

Tests::Tests(Basis *B) {
    b = B;
} // end constructor

Tests::~Tests() {
    delete b;
} // end deconstructor

void Tests::run_tests(int t) {
} // end function run_tests

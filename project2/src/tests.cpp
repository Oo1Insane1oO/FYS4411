//////////////////////////////////////////////////////////////////////////////
// Class containing functions for testinc class Basis.						//
//																			//
// All functions retunr true in case where the tests are same and false		//
// otherwise. Except for run_test which just runs all the functions.		//
//////////////////////////////////////////////////////////////////////////////

#include "tests.h" // header
#include <iostream>

Tests::Tests(Basis *B) {
    b = B;
} // end constructor

Tests::~Tests() {
    delete b;
} // end deconstructor

bool Tests::test_energies() {
    /* check that energy in state is correct */
    bool t = false;
    for (unsigned int i = 0; i < b->states.size(); ++i) {
        if (*(b->states)[i][4] == *(b->states)[i][0]+*(b->states)[i][1] + 1) {
            t = true;
        } else {
            t = false;
        } // emd ifelse
        if (!t) {
            break;
        } // end if
    } // end fori
    return t;
} // end function test_energies

// bool Tests::test_ratio() {
// } // end function test_ratio

void Tests::run_tests(int t) {
    /* run all tests and exit */
    if (t) {
        if(test_energies()) {
            std::cout << "Energies good" << std::endl;
        } else {
            std::cout << "Energies wrong" << std::endl;
        } // end ifelse
        if (t==2) {
            b->printStates();
        } // end ift
        exit(1);
    } // end if t
} // end function run_tests

//////////////////////////////////////////////////////////////////////////////
// Class containing functions for testinc class Basis.						//
//																			//
// All functions retunr true in case where the tests are same and false		//
// otherwise. Except for run_test which just runs all the functions.		//
//////////////////////////////////////////////////////////////////////////////

#include "tests.h" // header
#include <iostream>

Tests::Tests(Basis *B, VMC *V) {
    b = B;
    v = V;
    m = new Methods();
} // end constructor

Tests::~Tests() {
    delete b;
    delete v;
    delete m;
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

bool Tests::test_2particle() {
    /* check that energy in case of unperturbed harmonic oscillator system with
     * 2 electrons is correct */
    v->calculate(false);
    return ((m->variance(v->energy, v->energySq)) < 1e-14 ? true : false);
} // end function test_2particle

bool Tests::test_determinantratio() {
    /* check that ratio between determinants are correct in member function
     * determinantRatio of class methods */
    Eigen::Matrix3d A, B;
    A << 11, 2, 3,
         4, 5, 6,
         7, 78, 19;
    B << 11, 2, 3,
         14, 15, 16,
         7, 78, 19;
    return (std::fabs(m->determinantRatio(B,A.inverse(),1) -
                B.determinant()/A.determinant())<=1e10 ? true : false);
} // end function test_determinantratio

void Tests::run_tests(int t) {
    /* run all tests and exit */
    if (t) {
        if(test_energies()) {
            std::cout << "Energies good" << std::endl;
        } else {
            std::cout << "Energies wrong" << std::endl;
        } // end ifelse
        if(test_2particle()) {
            std::cout << "Energy unperturbed 2 electron good" << std::endl;
        } else { 
            std::cout << "Energy unperturbed 2 electron wrong" << std::endl;
        } // end ifelse
        if(test_determinantratio()) {
            std::cout << "Determinant ratio good" << std::endl;
        } else { 
            std::cout << "Determinant ratio wrong" << std::endl;
        } // end ifelse
        if (t==2) {
            b->printStates();
        } // end ift
        exit(1);
    } // end if t
} // end function run_tests

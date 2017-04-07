//////////////////////////////////////////////////////////////////////////////
// Class containing functions for testing class Methods and Basis.          //
//                                                                          //
// All functions return true in case where the tests are same and false     //
// otherwise. Except for run_test which just runs all the functions.        //
//////////////////////////////////////////////////////////////////////////////

#include "tests.h"
tests#include <iostream>

Tests::Tests(double w, int E, int n) {
    /* set basis and create Methods object */
    b = new Basis(w, E);
    meth = new Methods();

    // set number of particles
    A = n;
} // end constructor

Tests::~Tests() {
    delete b;
    delete meth;
} // end deconstructor

bool Tests::test_hermite() {
    /* test first 10 hermite polynomial for x=0 */
    double eps = 1e-12;
    bool t = false;
    std::vector<int> exact({1,0,-2,0,12,0,-120,0,1680,0,-30240});
    for (unsigned i = 0; i < exact.size(); ++i) {
        t = (std::fabs(meth->hermite(0,i)-exact[i])<=eps) ? true : false;
        if (!t) {
            break;
        } // end if
    } // end fori
    return (t ? 1 : 0);
} // end function test_hermite

bool Tests::test_convert() {
    /* check that energies are still same after converting from cartesian to
     * polar quantum numbers */
    bool t = false;
    std::vector<int> nm(2);
    for (unsigned int i = 0; i < b->n.size(); ++i) {
        for (unsigned int j = 0; j < b->n.size(); ++j) {
            nm = b->convertToPolar(b->n[i],b->n[j]);
            if((b->n[i] + b->n[j]) == (2*nm[0] + std::abs(nm[1]))) {
                t = true;
            } else {
                t = false;
                break;
            } // end ifelse
        } // end forj
        if (!t) {
            break;
        } // end if
    } // end fori
    return t;
} // end function test_convert

bool Tests::test_energies() {
    /* check that energy in state is correct */
    bool t = false;
    for (unsigned int i = 0; i < b->states.size(); ++i) {
        if (*(b->states)[i][4] == *(b->states)[i][0] + *(b->states)[i][1] + 1) {
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

bool Tests::test_hartreefock() {
    /* check that in case of no Coulomb interaction that function HartreeFock
     * reproduces the unperturbed energies. */
    int maxIterations = 1000;

    // set unperturbed interaction elements
    b->assemble(0,true);
    b->HartreeFock(A,maxIterations,1e-10);
    for (unsigned int i = 0; i < b->singleParticleEnergiesHartreeFock.size();
            ++i) {
        if(std::fabs(b->eps0Integrals[i]-b->singleParticleEnergiesHartreeFock[i])
                >= 1e-16) {
            return false;
        } // end if
    } // end fori
    return true;
} // end function test_hartreefock

void Tests::run_tests(int t) {
    /* run all tests and exit */
    if (t) {
        if(test_hermite()) {
            std::cout << "Hermite good" << std::endl;
        } else {
            std::cout << "Hermite wrong" << std::endl;
        } // end ifelse
        if(test_convert()) {
            std::cout << "Convert good" << std::endl;
        } else {
            std::cout << "Convert wrong" << std::endl;
        } // end ifelse
        if(test_energies()) {
            std::cout << "Energies good" << std::endl;
        } else {
            std::cout << "Energies wrong" << std::endl;
        } // end ifelse
        if(test_hartreefock()) {
            std::cout << "Hartree Fock unperturbed good" << std::endl;
        } else {
            std::cout << "Hartree Fock unperturbed good" << std::endl;
        } // end ifelse
        if (t==2) {
            b->printStates();
        } // end ift
        exit(1);
    } // endif
} // end function run_tests

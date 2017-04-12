#include "tests.h" // test functions
#include "vmc.h" // class basis
#include <stdlib.h> // atoi
#include <iostream> // cout
#include <chrono> // timer

//////////////////////////////////////////////////////////////////////////////
// Main file for running vmc algorithm                                      //
//////////////////////////////////////////////////////////////////////////////

int main(int argc, const char** argv) {
    /* main */

    if (argc < 5) {
        /* Print usage if number of command line arguments are to few */
        std::cout << 
            "USAGE: ./main 'omega' 'cutoff' 'particles' 'iterations' 'tests'" 
            << std::endl;
        std::cout <<
            "    " << "omega: (float) HO frequency\n" <<
            "    " << "cutoff: (num) Number of single particle orbitals\n" <<
            "    " << "particles: (num) Fermi level(closed shell)\n" <<
            "    " << "iterations: (num) Max iterations in VMC(MC cycles) \n" <<
            "    " << "tests: (1/0) indicating to run tests or not" <<
            std::endl;
        exit(1);
    } // end if

    // grab parameters as command line arguments
    double omega = atof(argv[1]);
    int Ec = atoi(argv[2]);
    int num = atoi(argv[3]);
    int maxIterations = atoi(argv[4]);
    int t = atoi(argv[5]);
    
    // set basis (cartesian)
    Basis *b = new Basis(omega, Ec);
    if (num>*(b->states[b->states.size()-1][5])) {
        std::cout << "Increase cutoff" << std::endl;
        exit(1);
    } // end if
    
    if (t) {
        /* run tests */
        Tests testObj = Tests(b);
        testObj.run_tests(t);
        exit(1);
    } // end if

    // set vmc object for calculations
    VMC vmcObj = VMC(b);
    Eigen::MatrixXd r = Eigen::MatrixXd::Zero(num,2);
    r(0,0) = 1;
    r(0,1) = -2;
    r(1,0) = 3;
    r(1,1) = 2;
    std::cout << b->trialWaveFunction(r,1.1,0.2,0) << std::endl;

    return 0;
} // end main

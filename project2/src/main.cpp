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
            "    " << "particles: (int) Fermi level(closed shell)\n" <<
            "    " << "iterations: (int) Max iterations in VMC(MC cycles) \n" <<
            "    " << "step: (float) step size in VMC \n" <<
            "    " << "tests: (1/0) indicating to run tests or not" <<
            std::endl;
        exit(1);
    } // end if

    // grab parameters as command line arguments
    double omega = atof(argv[1]);
    int num = atoi(argv[2]);
    int maxIterations = atoi(argv[3]);
    double step = atof(argv[4]);
    int t = atoi(argv[5]);

    if (num > 2) {
        std::cout << "Not implemented yet" << std::endl;
    } // end if
    
    // set basis (cartesian)
    Basis *b = new Basis(omega, num/2);
    
    if (t) {
        /* run tests */
        Tests testObj = Tests(b);
        testObj.run_tests(t);
        exit(1);
    } // end if

    // set vmc object for calculations
    VMC vmcObj = VMC(b,1,1);
    vmcObj.calculate(step, maxIterations);
    std::cout << "Energy: " << vmcObj.energy << std::endl;

    return 0;
} // end main

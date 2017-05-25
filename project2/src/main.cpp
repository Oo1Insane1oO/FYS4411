#include "tests.h" // test functions
#include "vmc.h" // class basis
#include <stdlib.h> // atoi
#include <iostream> // cout
#include <chrono> // timer
#include <iomanip> // setprecision
#include <algorithm> // find

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
            "    " << "importance sampling: (1/0) indicating to run with importance sampling or not" <<
            std::endl;
        exit(1);
    } // end if

    // grab parameters as command line arguments
    double omega = atof(argv[1]);
    int num = atoi(argv[2]);
    unsigned int maxIterations = atoi(argv[3]);
    double step = atof(argv[4]);
    int t = atoi(argv[5]);
    bool imp = atoi(argv[6]);
    
    // set basis (cartesian)
    Basis *b = new Basis(omega, num/2);

    // make sure number of particles is a magic number(closed shell)
    std::vector<int> magicNumber = b->getMagicNumbers();
    std::vector<int>::iterator it;
    it = std::find(magicNumber.begin(), magicNumber.end(), num);
    if (it == magicNumber.end()) {
        std::cout << "make sure num is a magic number N=2,6,12,20,30,42..." <<
            std::endl;
        exit(1);
    } // end if

    std::cout << "Basis made" << std::endl;
    
    // set vmc object for calculations
//     VMC *vmcObj = new VMC(b,1,0.4,2,step,maxIterations,imp);
    VMC *vmcObj = new VMC(b,1,0.4,2,step,maxIterations,imp);
//     VMC *vmcObj = new VMC(b,1.03741,0.472513,2,step,maxIterations,imp);
    
    if (t) {
        /* run tests */
        Tests testObj = Tests(b,vmcObj,num);
        testObj.run_tests(t);
        exit(1);
    } // end if

    // run calculations
    vmcObj->calculate(imp);
    std::cout << std::setprecision(10) << "<E> = " << vmcObj->energy << ", " <<
        "<E^2> = " << vmcObj->energySq << std::endl;
    std::cout << std::setprecision(10) << "<E^2> - <E>^2 = " <<
        (vmcObj->energySq - pow(vmcObj->energy,2)) << std::endl;

    std::cout << "alpha: " << vmcObj->alpha << ", beta: " << vmcObj->beta <<
        std::endl;

    // free objects
    delete b;
    delete vmcObj;

    return 0;
} // end main

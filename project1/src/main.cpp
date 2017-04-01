#include "tests.h" // test functions
#include "basis.h" // class basis
#include <stdlib.h> // atoi
#include <iostream> // cout
#include <chrono> // timer

//////////////////////////////////////////////////////////////////////////////
// Main file for running Hartree-Fock algorithm                             //
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
            "    " << "iterations: (num) Max iterations in Hartree-Fock\n" <<
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

    if (t) {
        /* run tests */
        Tests testObj = Tests(omega, Ec, num);
        testObj.run_tests(t);
    } // end if

    // set basis (cartesian)
    Basis b = Basis(omega, Ec);
    if (num>b.M[b.M.size()-1]) {
        std::cout << "Increase Ec" << std::endl;
        exit(1);
    } // end if

    // assemble integrals
    std::chrono::steady_clock::time_point begin;
    begin = std::chrono::steady_clock::now();
    b.assemble();
    std::chrono::steady_clock::time_point end;
    end = std::chrono::steady_clock::now();
    std::cout << "Assemble time: " <<
        std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() <<
        std::endl;

    // run Hartree-Fock algorithm
    begin = std::chrono::steady_clock::now();
    b.HartreeFock(num,maxIterations,1e-10);
    end = std::chrono::steady_clock::now();
    std::cout << std::fixed << "Single particles energies:" << std::endl;
//     for (unsigned int i = 0; i < b.M.size(); ++i) {
//         std::cout << b.singleParticleEnergiesHartreeFock[b.M[i]-1] << " " <<
//             std::endl;
//         if (b.M[i]==static_cast<int>(num)) {
//             std::cout << b.singleParticleEnergiesHartreeFock[num] <<
//                 std::endl;
//         } // end if
//     } // end fori
    for (unsigned int i = 0; i < b.singleParticleEnergiesHartreeFock.size();
            ++i) {
        std::cout << i+1 << " " << b.singleParticleEnergiesHartreeFock[i] <<
            std::endl;
    } // end fori
    std::cout << "Hartree-Fock energy after " << maxIterations << 
        " Iterations: " << b.E0HartreeFock << std::endl;
    std::cout << "HF time: " <<
        std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() <<
        std::endl;
    return 0;
} // end main

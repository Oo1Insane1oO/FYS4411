//////////////////////////////////////////////////////////////////////////////
// Class for calculating variational Monte Carlo energies.                  //
//                                                                          //
// Main Functions:                                                          //
// See the individual functions for specific behavior.                      //
//////////////////////////////////////////////////////////////////////////////

#include "vmc.h" // header
#include <iostream>
#include <fstream>
#include <math.h>
#include <random>

VMC::VMC(Basis *B) {
    b = B;
    meth = new Methods(); 
} // end constructor

VMC::~VMC() {
    delete meth;
    delete b;
} // end deconstructor
        
double VMC::localEnergy2(Eigen::MatrixXd r1, Eigen::MatrixXd r2, bool coulomb) {
    /* calculate analytic expression of local energy for 2 electrons */
    double r12 = (r1-r2).norm();
    double denom = 1 + beta*r12;
    return 0.5 * pow(b->omega,2) * (pow(alpha,2) + 1) * (r1.squaredNorm() +
            r2.squaredNorm()) - 2*alpha*b->omega + a/pow(denom,2) *
        ((a/pow(denom,2) - alpha*b->omega*r12 + 1/r12 - 2*beta/denom)) +
        (coulomb ? 1/r12 : 0);
} // end function localEnergy

void calculate() {
} // end function calculate

double VMC::metropolisTest(double densityRatio, double proposedRatio) {
    /* perform metropolist test */
    std::mt19937_64 randomGenerator();
    return meth->min(1.,densityRatio*proposedRatio);
} // end function metropolisTest

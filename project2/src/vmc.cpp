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

VMC::VMC(Basis *B, double alp, double bet) {
    alpha = alp;
    beta = bet;
    a = 1;
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

void VMC::initialize(unsigned long int seed) {
    /* initialize positions R */
    std::mt19937_64 mt(seed); 
    std::uniform_real_distribution<double> randomReal(-1,1);
    R.resize(b->ECut,b->ECut);
    for (int i = 0; i < b->ECut; ++i) {
        for (int j = 0; j < 2; ++j) {
            R(i,j) = randomReal(mt);
        } // end forj
    } // end fori
} // end initialize positions

void VMC::calculate(double step, int cycles) {
    a = 0; //TODO: FIX THIS
    std::mt19937_64 mt(85456);
    std::uniform_real_distribution<double> r(0,1);
    Eigen::MatrixXd Rp;
    initialize();
    double P, Pp;
    for (int i = 0; i < cycles; ++i) {
        P = b->trialWaveFunction(R,alpha,beta,a);
        Rp = (R.array() + r(mt) * step).matrix();
        Pp = b->trialWaveFunction(Rp,alpha,beta,a);
        if (metropolisTest(Pp/P,1)>=1) {
            R = Rp;
        } // end if
    } // end fori
} // end function calculate

double VMC::metropolisTest(double densityRatio, double proposedRatio) {
    /* perform metropolist test */
    return meth->min(1.,densityRatio*proposedRatio);
} // end function metropolisTest

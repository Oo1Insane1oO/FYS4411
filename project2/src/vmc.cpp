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
        
double VMC::localEnergy2(Eigen::MatrixXd R, bool coulomb) {
    /* calculate analytic expression of local energy for 2 electrons */
    Eigen::MatrixXd r1 = R.block<1,2>(0,0);
    Eigen::MatrixXd r2 = R.block<1,2>(1,0);
    double r12 = (r1-r2).norm();
    double denom = 1 + beta*r12;
    return 0.5 * pow(b->omega,2) * (1 - pow(alpha,2)) * (r1.squaredNorm() +
            r2.squaredNorm()) + 2*alpha*b->omega - a/pow(denom,2) *
        ((a/pow(denom,2) + alpha*b->omega*r12 + 1/r12 - 2*beta/denom)) +
        (coulomb ? 1/r12 : 0);
} // end function localEnergy

void VMC::initialize(unsigned long int seed) {
    /* initialize positions R */
    a = 0; //TODO: FIX THIS
    std::mt19937_64 mt(seed); 
    std::uniform_real_distribution<double> randomReal(-1,1);
    R.resize(2*b->ECut,2*b->ECut);
    for (int i = 0; i < b->ECut; ++i) {
        for (int j = 0; j < 2; ++j) {
            R(i,j) = randomReal(mt);
        } // end forj
    } // end fori
    energy = 0;
} // end initialize positions

void VMC::calculate(double step, int cycles) {
    std::mt19937_64 mt(85456);
    std::uniform_real_distribution<double> r(0,1);
    Eigen::MatrixXd Rp;
    initialize();
    double P;
    for (int i = 0; i < cycles; ++i) {
        P = pow(b->trialWaveFunction(R,alpha,beta,a),2);
        energy += P*localEnergy2(R,false);
        Rp = (R.array() + r(mt) * step).matrix();
        if (metropolisTest(pow(b->trialWaveFunction(Rp,alpha,beta,a),2)/P,1)) {
            R = Rp;
        } // end if
    } // end fori
    energy /= cycles;
} // end function calculate

bool VMC::metropolisTest(double densityRatio, double proposedRatio) {
    /* perform metropolist test */
    return (meth->min(1.,densityRatio*proposedRatio)>=1 ? true : false);
} // end function metropolisTest

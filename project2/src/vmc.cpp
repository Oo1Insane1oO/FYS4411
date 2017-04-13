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
    double denomsq = denom*denom;
    return 0.5 * pow(b->omega,2) * (1 - pow(alpha,2)) * (r1.squaredNorm() +
            r2.squaredNorm()) + 2*alpha*b->omega - a/denomsq *
        ((a/denomsq + alpha*b->omega*r12 + 1/r12 - 2*beta/denom)) +
        (coulomb ? 1/r12 : 0);
} // end function localEnergy

void VMC::initialize(unsigned long int seed, double bound) {
    /* initialize positions R */
    a = 0; //TODO: FIX THIS
    std::mt19937_64 mt(seed); 
    std::uniform_real_distribution<double> randomReal(-bound,bound);
    R.resize(2*b->ECut,2*b->ECut);
    for (int i = 0; i < b->ECut; ++i) {
        for (int j = 0; j < 2; ++j) {
            R(i,j) = randomReal(mt);
        } // end forj
    } // end fori
    energy = 0;
} // end initialize positions

void VMC::calculate(double step, int maxIterations, unsigned long int seed) {
    std::mt19937_64 mt(seed);
    std::uniform_real_distribution<double> dist(0,1);
    initialize();
    Eigen::MatrixXd Rp = Eigen::MatrixXd::Zero(R.rows(),R.cols());
    double P, r;
    int cycles = 0;
    while (cycles < maxIterations) {
        P = pow(b->trialWaveFunction(R,alpha,beta,a),2);
        r = dist(mt);
        Rp = (R.array() + r * step).matrix();
        if (pow(b->trialWaveFunction(Rp,alpha,beta,a),2) / P >= r) {
            R = Rp;
        } // end if
        energy += localEnergy2(R,false);
        cycles++;
    } // end while
    energy /= cycles;
} // end function calculate

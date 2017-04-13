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

VMC::VMC(Basis *B, double alp, double bet, unsigned int d) {
    alpha = alp;
    beta = bet;
    a = 1;
    b = B;
    dim = d;
    meth = new Methods(); 
} // end constructor

VMC::~VMC() {
    delete meth;
    delete b;
} // end deconstructor
        
double VMC::localEnergy2(Eigen::MatrixXd R, bool coulomb) {
    /* calculate analytic expression of local energy for 2 electrons */
    double r12 = (R.row(0) - R.row(1)).norm();
    double denom = 1 + beta*r12;
    double denomsq = denom*denom;
    return 0.5 * pow(b->omega,2) * (1 - pow(alpha,2)) * (R.row(0).squaredNorm()
            + R.row(1).squaredNorm()) + 2*alpha*b->omega - a/denomsq *
        ((a/denomsq + alpha*b->omega*r12 + 1/r12 - 2*beta/denom)) + (coulomb ?
            1/r12 : 0);
} // end function localEnergy

void VMC::calculate(double step, int maxIterations, unsigned long int seed) {
    /* function for running Monte Carlo integration */

    // initialize Mersenne Twister random number generator and uniform
    // distribution engine
    std::mt19937_64 mt(seed);
    std::uniform_real_distribution<double> dist(0,1);

    // initialize position
    Eigen::MatrixXd oldPositions = Eigen::MatrixXd::Zero(2*b->ECut, dim);
    Eigen::MatrixXd newPositions = Eigen::MatrixXd::Zero(2*b->ECut, dim);
    for (unsigned int i = 0; i < oldPositions.rows(); ++i) {
        for (unsigned int j = 0; j < oldPositions.cols(); ++j) {
            oldPositions(i,j) = step * (dist(mt)-0.5);
        } // end forj
    } // end fori
    energy = 0;
    energySq = 0;
    newPositions = oldPositions;

    int cycles = 0;
    double Pnew, Pold, tmpEnergy;
    while (cycles < maxIterations) {
        /* run Monte Carlo cycles */
        // set current wave function
        Pold = pow(b->trialWaveFunction(oldPositions,alpha,beta,a),2);
        for (unsigned int i = 0; i < oldPositions.rows(); ++i) {
            /* loop over number of particles */
            for (unsigned int j = 0; j < dim; ++j) {
                /* propose new position */
                newPositions(i,j) = oldPositions(i,j) + step*(dist(mt)-0.5);
            } // end forj

            // calculate new PDF
            Pnew = pow(b->trialWaveFunction(newPositions,alpha,beta,a),2);

            if (Pnew/Pold >= dist(mt)) {
                /* update positions according to Metropolis test */
                oldPositions.row(i) = newPositions.row(i);
                Pold = Pnew;
            } else {
                /* reset position */
                newPositions.row(i) = oldPositions.row(i);
            } // end if

            // update energy and increment cycles
            tmpEnergy = localEnergy2(newPositions,false);
            energy += tmpEnergy;
            energySq += tmpEnergy*tmpEnergy;
        } // end fori
        cycles++;
    } // end while

    // calculate final energy estimation
    energy /= cycles * newPositions.rows();
    energySq /= cycles * newPositions.rows();
} // end function calculate

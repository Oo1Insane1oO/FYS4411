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

VMC::VMC(Basis *B, double alp, double bet, unsigned int d, double s, unsigned
        int max, bool sample) {
    alpha = alp;
    beta = bet;
    a = 0;
    b = B;
    dim = d;
    step = s;
    dt = 0.001;
    maxIterations = max;
    imp = sample;
    meth = new Methods(); 
} // end constructor

VMC::~VMC() {
    delete meth;
} // end deconstructor
        
double VMC::localEnergy2(Eigen::MatrixXd R, bool coulomb) {
    /* calculate analytic expression of local energy for 2 electrons */
    double r12 = (R.row(0) - R.row(1)).norm();
    double denom = 1 + beta*r12;
    double denomsq = denom*denom;
    return 0.5 * pow(b->omega,2) * (1 - pow(alpha,2)) * (R.row(0).squaredNorm()
            + R.row(1).squaredNorm()) + 2*alpha*b->omega - a/denomsq *
        ((a/denomsq - alpha*b->omega*r12 + 1/r12 - 2*beta/denom)) + (coulomb ?
            1/r12 : 0);
} // end function localEnergy

void VMC::diff(Eigen::MatrixXd R, Eigen::MatrixXd &der) {
    /* calculate first derivative ratio of single particle wave functions */
    double r12 = (R.row(0) - R.row(1)).norm();
    double denom = 1 + beta*r12;
    for (unsigned int i = 0; i < R.rows(); ++i) {
        for (unsigned int j = 0; j < R.cols(); ++j) {
            der(i,j) = -alpha*b->omega*R(i,j) + a*(R(i,j) -
                    R(j,i))/(r12*denom*denom);
        } // end forj
    } // end fori
} // end function

// double VMC::localEnergyDiff(Eigen::MatrixXd R, bool coulomb) {
//     /* calculate analytic expression of local energy for 2 electrons */
//     double dx = 0.0001;
//     return 0.5 * (-diff2(R,dx) + pow(b->omega,2) * (R.row(0).squaredNorm() +
//                 R.row(1).squaredNorm())) + (coulomb ?
//             1/(R.row(0)-R.row(1)).norm() : 0);
// } // end function localEnergyDiff

// double VMC::diff2(Eigen::MatrixXd R, double dx) {
//     /* calculate second derivative for all positions in x using central
//      * difference scheme */
//     double diff = 0;
//     double tmpDiff;
//     Eigen::MatrixXd Rpm = R;
//     double mid = 2*b->trialWaveFunction(R,alpha,beta,a);
//     for (unsigned int i = 0; i < R.rows(); ++i) {
//         for (unsigned int j = 0; j < R.cols(); ++j) {
//             Rpm(i,j) += dx;
//             tmpDiff = b->trialWaveFunction(Rpm,alpha,beta,a) - mid;
//             Rpm(i,j) -= 2*dx;
//             tmpDiff += b->trialWaveFunction(Rpm,alpha,beta,a);
//             diff += tmpDiff / (dx*dx);
//             Rpm = R;
//         } // end forj
//     } // end fori
//     return diff;
// } // end function diff2

void VMC::setSeed(unsigned long int s) {
    /* set seed */
    seed = s;
} // end function setSeed

unsigned long int VMC::getSeed() {
    /* get seed */
    return seed;
} // end function setSeed

void VMC::calculate(bool unperturb) {
    /* function for running Monte Carlo integration */

    // initialize Mersenne Twister random number generator and uniform
    // distribution engine
    std::mt19937_64 mt(seed);
    std::uniform_real_distribution<double> dist(0,1);
    std::normal_distribution<double> normDist(1,0);

    // initialize position
    Eigen::MatrixXd oldPositions = Eigen::MatrixXd::Zero(2*b->ECut, dim);
    Eigen::MatrixXd newPositions = Eigen::MatrixXd::Zero(2*b->ECut, dim);
    for (unsigned int i = 0; i < oldPositions.rows(); ++i) {
        for (unsigned int j = 0; j < oldPositions.cols(); ++j) {
            if (imp) {
                oldPositions(i,j) = normDist(mt) * sqrt(dt);
            } else {
                oldPositions(i,j) = step * (dist(mt)-0.5);
            } // end ifelse
        } // end forj
    } // end fori
    energy = 0;
    energySq = 0;
    newPositions = oldPositions;

    Eigen::MatrixXd qForceOld, qForceNew;
    if (imp) {
        qForceOld = Eigen::MatrixXd::Zero(oldPositions.rows(),oldPositions.cols());
        qForceNew = Eigen::MatrixXd::Zero(oldPositions.rows(),oldPositions.cols());
    } // end if

    double testRatio;
    unsigned int cycles = 0;
    double tmpEnergy, greensFunctionRatio;
    Eigen::MatrixXd oldWaveFunction, newWaveFunction, oldInverse;
    Eigen::MatrixXd newInverse = Eigen::MatrixXd::Zero(oldPositions.rows(),
            oldPositions.rows());
    while (cycles < maxIterations) {
        /* run Monte Carlo cycles */
        oldWaveFunction = b->trialWaveFunction(oldPositions,alpha,beta,a);
        oldInverse = oldWaveFunction.inverse();
        if (imp) {
            diff(oldPositions,qForceOld);
        } // end if
        for (unsigned int i = 0; i < oldPositions.rows(); ++i) {
            /* loop over number of particles */
            for (unsigned int j = 0; j < dim; ++j) {
                /* propose new position */
                if (imp) {
                    newPositions(i,j) = oldPositions(i,j) +
                        0.5*qForceOld(i,j)*dt + normDist(mt)*sqrt(dt);
                } else {
                    newPositions(i,j) = oldPositions(i,j) + step*(dist(mt)-0.5);
                } // end ifelse
            } // end forj

            // calculate new PDF (probability distribution function)
            newWaveFunction = b->trialWaveFunction(newPositions,alpha,beta,a);
            if (imp) {
                diff(newPositions, qForceNew);
            } // end if

            // calculate Greens funcition ratio
            if (imp) {
                greensFunctionRatio = 0;
                for (unsigned int j = 0; j < oldPositions.rows(); ++j) {
                    greensFunctionRatio -= pow(oldPositions(i,j) -
                            normDist(mt)*sqrt(dt),2) - pow(newPositions(i,j) -
                            normDist(mt)*sqrt(dt),2);
                } // end forj
                greensFunctionRatio = exp(greensFunctionRatio);
            } // end if

            testRatio = pow(meth->determinantRatio(newWaveFunction, oldInverse, i), 2);
            if (imp) {
                testRatio *= greensFunctionRatio;
            } //end if

            if (testRatio >= dist(mt)) {
                /* update positions according to Metropolis test */
                oldPositions.row(i) = newPositions.row(i);
                oldWaveFunction = newWaveFunction;
                qForceOld = qForceNew;
            } else {
                /* reset position */
                newPositions.row(i) = oldPositions.row(i);
                qForceNew = qForceOld;
            } // end if

            // update energy and increment cycles
            tmpEnergy = localEnergy2(newPositions,unperturb);
//             tmpEnergy = localEnergyDiff(newPositions,false);
            energy += tmpEnergy;
            energySq += tmpEnergy*tmpEnergy;
            meth->updateMatrixInverse(oldWaveFunction, newWaveFunction,
                    oldInverse, newInverse, i);
            oldInverse = newInverse;
        } // end fori
        cycles++;
    } // end while

    // calculate final energy estimation
    energy /= cycles * newPositions.rows();
    energySq /= cycles * newPositions.rows();
} // end function calculate

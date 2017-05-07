//////////////////////////////////////////////////////////////////////////////
// Class for calculating variational Monte Carlo energies.                  //
//                                                                          //
// Main Functions:                                                          //
// See the individual functions for specific behavior.                      //
//////////////////////////////////////////////////////////////////////////////

#include "vmc.h" // header
#include "hermite.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <random>

VMC::VMC(Basis *B, double alp, double bet, unsigned int d, double s, unsigned
        int max, bool sample) {
    alpha = alp;
    beta = bet;
    b = B;
    dim = d;
    step = s;
    maxIterations = max;
    imp = sample;

    meth = new Methods(); 
} // end constructor

VMC::~VMC() {
    delete meth;
} // end deconstructor
        
// double VMC::localEnergy2(const Eigen::MatrixXd &R, bool coulomb) {
//     /* calculate analytic expression of local energy for 2 electrons */
//     double r12 = (R.row(0) - R.row(1)).norm();
//     double denom = 1 + beta*r12;
//     double denomsq = denom*denom;
//     return 0.5 * pow(b->omega,2) * (1 - alpha*alpha) * (R.row(0).squaredNorm()
//             + R.row(1).squaredNorm()) + 2*alpha*b->omega + (coulomb ?
//             -1/denomsq * (1/denomsq - alpha*b->omega*r12 + 1/r12 -
//                     2*beta/denom) + 1/r12 : 0);
// } // end function localEnergy

double VMC::localEnergy2(const Eigen::MatrixXd &R, bool coulomb) {
    /* calculate analytic expression of local energy */
    double nx, ny, nxHermiteFactor, nyHermiteFactor, rk, rkj, jFactor, denom,
           a;
    double E = 0;
    for (unsigned int k = 0; k < R.rows(); ++k) {
        /* loop over particles */
        rk = R.row(k).norm();
        nx = *(b->states[k][0]);
        ny = *(b->states[k][1]);
        nxHermiteFactor = nx*(nx-1)*H(R(k,0),nx-2)/H(R(k,0),nx);
        nyHermiteFactor = ny*(ny-1)*H(R(k,1),ny-2)/H(R(k,1),ny);
        E += 0.5 * pow(b->omega,2)*rk*rk;
        E -= b->omega*(2-alpha) * (nxHermiteFactor + nyHermiteFactor) +
            0.5*alpha*b->omega * (alpha*b->omega*rk*rk - 2*(nx+ny+1));
        if (coulomb) {
            /* Add Jastrow part */
            for (unsigned int j = 0; j < R.rows(); ++j) {
                if (j != k) {
                    a = b->padejastrow(k,j);
                    rkj = (R.row(k) - R.row(j)).norm();
                    denom = 1 + beta*rkj;
                    jFactor = 0.5*a/pow(denom,2); 
                    E -= jFactor * (2/rkj*(((nx + nxHermiteFactor)/R(k,0) -
                                    alpha*b->omega*R(k,0))*(R(k,0)-R(j,0)) +
                                ((ny + nyHermiteFactor)/R(k,1) -
                                 alpha*b->omega*R(k,1))*(R(k,1)-R(j,1))) +
                            1/rkj - 2*beta/denom + a / pow(denom,2));
                    if (j > k) {
                        /* Coulomb part */
                        E += 1/rkj;
                    } // end if
                } // end if
            } // end forj
        } // end if
    } // end fork
    return E;
} // end function localEnergy

// void VMC::diff(const Eigen::MatrixXd &R, Eigen::MatrixXd &der) {
//     /* calculate first derivative ratio of single particle wave functions */
//     double r12 = (R.row(0) - R.row(1)).norm();
//     double denom = r12 * pow(1+beta*r12, 2);
//     for (unsigned int i = 0; i < R.rows(); ++i) {
//         for (unsigned int j = 0; j < R.cols(); ++j) {
//             der(i,j) = -alpha*b->omega*R(i,j) + b->padejastrow(i,j)*(R(i,j) -
//                     R(i+((i%2 || i==1) ? -1 : 1),j))/denom;
//         } // end forj
//     } // end fori
// } // end function

void VMC::diff(const Eigen::MatrixXd &R, Eigen::MatrixXd &der) {
    /* calculate first derivative ratio of single wave functions */
    double nx, ny, nxHermiteFactor, nyHermiteFactor, tmpkVal;
    for (unsigned int k = 0; k < R.rows(); ++k) {
        nx = *(b->states[k][0]);
        nx = *(b->states[k][1]);
        nxHermiteFactor = nx*(nx-1)*H(R(k,0),nx-2)/H(R(k,0),nx);
        nyHermiteFactor = ny*(ny-1)*H(R(k,1),ny-2)/H(R(k,1),ny);
        tmpkVal = (nx + nxHermiteFactor) / R(k,0) + (ny + nyHermiteFactor) /
            R(k,1) - alpha*b->omega * (R(k,0) + R(k,1));
        for (unsigned int j = 0; j < R.rows(); ++j) {
            der(k,j) = tmpkVal;
            if (j != k) {
                der(k,j) += b->padejastrow(k,j) * (R(k,0) - R(j,0) + R(k,1) -
                        R(j,1)) / pow((1 + beta*(R.row(k)-R.row(j)).norm()),2);
            } // end if
        } // end forj
    } // end fork
} // end function diff

double VMC::localEnergyDiff(Eigen::MatrixXd &psiD, Eigen::MatrixXd &psiU, const
        Eigen::MatrixXd &R) {
    /* calculate local energy electrons */
    // set kinetic part and calculate potential
    double dx = 0.00005;
    double diff = -diff2(psiD, psiU, R, dx);
    for (unsigned int i = 0; i < R.cols(); ++i) {
        /* calculate potential part */
        diff += pow(b->omega,2) * R.row(i).squaredNorm();
    } // end fori
    return 0.5 * diff;
} // end function localEnergyDiff

double VMC::diff2(Eigen::MatrixXd &psiD, Eigen::MatrixXd &psiU, const
        Eigen::MatrixXd &R, double dx) {
    /* calculate second derivative for all positions in R using central
     * difference scheme */
    double diff = 0;
    double tmpDiff;
    Eigen::MatrixXd Rpm = R;
    b->setTrialWaveFunction(psiD,psiU,R,alpha);
    double mid = 2*psiD.determinant()*psiU.determinant();
    for (unsigned int i = 0; i < R.rows(); ++i) {
        for (unsigned int j = 0; j < R.cols(); ++j) {
            Rpm(i,j) += dx;
            b->setTrialWaveFunction(psiD,psiU,Rpm,alpha);
            tmpDiff = psiD.determinant() * psiU.determinant() - mid;
            Rpm(i,j) -= 2*dx;
            b->setTrialWaveFunction(psiD,psiU,Rpm,alpha);
            tmpDiff += psiD.determinant() * psiU.determinant();
            diff += tmpDiff / (dx*dx);
            Rpm = R;
        } // end forj
    } // end fori
    return diff;
} // end function diff2

void VMC::setSeed(unsigned long int s) {
    /* set seed */
    seed = s;
} // end function setSeed

unsigned long int VMC::getSeed() {
    /* get seed */
    return seed;
} // end function setSeed

void VMC::calculate(bool perturb) {
    /* function for running Monte Carlo integration */

    // initialize Mersenne Twister random number generator and uniform
    // distribution engine
    std::mt19937_64 mt(seed);
    std::uniform_real_distribution<double> dist(0,1);
    std::normal_distribution<double> normDist(0,1);

    // initialize position
    Eigen::MatrixXd oldPositions = Eigen::MatrixXd::Zero(2*b->ECut, dim);
    Eigen::MatrixXd newPositions = Eigen::MatrixXd::Zero(2*b->ECut, dim);
    for (unsigned int i = 0; i < oldPositions.rows(); ++i) {
        for (unsigned int j = 0; j < oldPositions.cols(); ++j) {
            if (imp) {
                oldPositions(i,j) = normDist(mt) * sqrt(step);
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
        qForceOld = Eigen::MatrixXd::Zero(oldPositions.rows(),
                oldPositions.cols());
        qForceNew = Eigen::MatrixXd::Zero(oldPositions.rows(),
                oldPositions.cols());
    } // end if

    unsigned int halfSize = oldPositions.rows()/2;
    double testRatio;
    double determinantRatioD = 0;
    double determinantRatioU = 0;
    unsigned int cycles = 0;
    double tmpEnergy, greensFunctionRatio;
    Eigen::MatrixXd oldD = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd oldU = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd newD = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd newU = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd oldInvD = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd oldInvU = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd newInvD = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd newInvU = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    while (cycles < maxIterations) {
        /* run Monte Carlo cycles */
        b->setTrialWaveFunction(oldD, oldU, oldPositions, alpha);
        oldInvD = oldD.inverse();
        oldInvU = oldU.inverse();
        if (imp) {
            diff(oldPositions,qForceOld);
            qForceOld *= 2;
        } // end if
        for (unsigned int i = 0; i < oldPositions.rows(); ++i) {
            /* loop over number of particles */
            for (unsigned int j = 0; j < dim; ++j) {
                /* propose new position */
                if (imp) {
                    newPositions(i,j) = oldPositions(i,j) +
                        0.5*qForceOld(i,j)*step + normDist(mt)*sqrt(step);
                } else {
                    newPositions(i,j) = oldPositions(i,j) +
                        step*(dist(mt)-0.5);
                } // end ifelse
            } // end forj

            // calculate new PDF (probability distribution function)
            b->setTrialWaveFunction(newD, newU, newPositions, alpha);
            if (imp) {
                /* set new quantum force */
                diff(newPositions, qForceNew);
                qForceNew *= 2;
            } // end if

            // calculate Greens function ratio
            if (imp) {
                greensFunctionRatio = exp(-(0.5*(qForceOld.row(i).array() +
                                qForceNew.row(i).array()) * (0.25*step *
                                (qForceOld.row(i).array() -
                                 qForceNew.row(i).array()) -
                                newPositions.row(i).array() +
                                oldPositions.row(i).array())).matrix().sum());
// 
//                 greensFunctionRatio = exp(((oldPositions.row(i).array() -
//                                 newPositions.row(i).array() - 0.5 * step *
//                                 qForceNew.row(i).array()).matrix().squaredNorm()
//                             - (newPositions.row(i).array() -
//                                 oldPositions.row(i).array() - 0.5
//                              * step *
//                              qForceOld.row(i).array()).matrix().squaredNorm())
//                         / (2*step));
            } // end if

            if ((i<halfSize) || (determinantRatioD==0)) {
                determinantRatioD = meth->determinantRatio(newD, oldInvD, i/2);
            } else if ((i>=halfSize) || (determinantRatioU==0)) {
                determinantRatioU = meth->determinantRatio(newU, oldInvU, i/2);
            } // end ifelseif

            testRatio = pow(determinantRatioD * determinantRatioU,2) *
                (!perturb ?  1 : exp(2*(b->jastrow(newPositions,beta) -
                                     b->jastrow(oldPositions,beta))));
            if (imp) {
                /* importance sampling */
                testRatio *= greensFunctionRatio;
            } //end if

            if (testRatio >= dist(mt)) {
                /* update positions according to Metropolis test */
                oldPositions.row(i) = newPositions.row(i);
                oldD = newD;
                oldU = newU;
                if (imp) {
                    qForceOld = qForceNew;
                } // end if
            } else {
                /* reset position */
                newPositions.row(i) = oldPositions.row(i);
                if (imp) {
                    qForceNew = qForceOld;
                } // end if
            } // end if

            // update energy and increment cycles
            tmpEnergy = localEnergy2(newPositions,perturb);
//             tmpEnergy = localEnergyDiff(newD,newU,newPositions) /
//                 (newD.determinant()*newU.determinant());
            energy += tmpEnergy;
            energySq += tmpEnergy*tmpEnergy;
            if (i < halfSize) {
                meth->updateMatrixInverse(oldD, newD, oldInvD, newInvD,
                        determinantRatioD, i/2);
            } else {
                meth->updateMatrixInverse(oldU, newU, oldInvU, newInvU,
                        determinantRatioU, i/2);
            } // end ifelse
        } // end fori
        cycles++;
    } // end while

    // calculate final energy estimation
    energy /= cycles * newPositions.rows();
    energySq /= cycles * newPositions.rows();
} // end function calculate

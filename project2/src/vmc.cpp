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
// 
// void VMC::diff(const Eigen::MatrixXd &R, Eigen::MatrixXd &der) {
//     /* calculate first derivative ratio of single particle wave functions */
//     double r12 = (R.row(0) - R.row(1)).norm();
//     double denom = r12 * pow(1+beta*r12, 2);
//     for (unsigned int i = 0; i < R.rows(); ++i) {
//         for (unsigned int j = 0; j < R.cols(); ++j) {
//             der(i,j) = -alpha*b->omega*R(i,j) + (R(i,j)-R(i+((i%2 || i==1) ?
//                             -1 : 1),j))/denom;
//         } // end forj
//     } // end fori
// } // end function
//  
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

void VMC::diff(const Eigen::MatrixXd &R, Eigen::MatrixXd &der) {
    /* calculate first derivative ratio of wave functions */
    for (unsigned int k = 0; k < R.rows(); ++k) {
        updateDiff(R,der,k);
    } // end fork
} // end function diff

void VMC::updateDiff(const Eigen::MatrixXd &R, Eigen::MatrixXd &der, unsigned
        int k) {
    /* calculate first derivative ratio of wave functions for partikle k */
    double rkj;
    for (unsigned int d = 0; d < R.cols(); ++d) {
        der(k,d) = (*(b->states[k][d])*(1 + (*(b->states[k][d])-1) *
                    H(R(k,d),*(b->states[k][d])-2) /
                    H(R(k,d),*(b->states[k][d]))))/R(k,d) -
            alpha*b->omega*R(k,d);
        for (unsigned int j = 0; j < R.rows(); ++j) {
            if (j != k) {
                rkj = (R.row(k) - R.row(j)).norm();
                der(k,d) += b->padejastrow(k,j) * (R(k,d)-R(j,d)) /
                    (rkj*pow(1+beta*rkj,2));
            } // end if
        } // end forj
    } // end ford
} // end function updateDiff

double VMC::localEnergyDiff(Eigen::MatrixXd &psiD, Eigen::MatrixXd &psiU, const
        Eigen::MatrixXd &R, bool coulomb) {
    /* calculate local energy electrons */
    double dx = 1e-5;
    
    // set kinetic part and calculate potential
    double diff = -diff2(psiD, psiU, R, dx) /
        (psiD.determinant()*psiU.determinant());
    for (unsigned int i = 0; i < R.rows(); ++i) {
        /* calculate potential part */
        diff += pow(b->omega,2) * R.row(i).squaredNorm();
    } // end fori
    diff *= 0.5;

    // calculate coulomb part
    if (coulomb) {
        for (unsigned int i = 0; i < R.rows(); ++i) {
            for (unsigned int j = i+1; j < R.rows(); ++j) {
                diff += 1/(R.row(i) - R.row(j)).norm();
            } // end forj
        } // end fori
    } // end if
    return diff;
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
            Rpm(i,j) += dx;
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
    double testRatio, tmpEnergy, tmpRsum, tmphsum, tmphsum3, tmpELalpRsum,
           tmpELbethsum, transitionRatio, denom, denomsq, denomcu, a, rkl,
           Hxfactor, Hyfactor;
    double ELhsum3 = 0;
    double hsum3 = 0;
    double Rsum = 0;
    double hsum = 0;
    double hsqsum = 0;
    double Rsqsum = 0;
    double ELR = 0;
    double ELh = 0;
    double ELRsq = 0;
    double ELhsq = 0;
    double ELRh = 0;
    double ELalpRsum = 0;
    double ELbethsum = 0;
    double determinantRatioD = 1;
    double determinantRatioU = 1;
    unsigned int cycles = 0;
    int nkx, nky;
    Eigen::MatrixXd oldD = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd oldU = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd newD = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd newU = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd oldInvD = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd oldInvU = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd newInvD = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd newInvU = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd HessenMatrix = Eigen::MatrixXd::Zero(2,2);
    Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(2,1);
    Eigen::MatrixXd newAlphaBeta = Eigen::MatrixXd::Zero(2,1);
    newAlphaBeta(0) = alpha;
    newAlphaBeta(1) = beta;
    b->setTrialWaveFunction(oldD, oldU, oldPositions, alpha);
    oldInvD = oldD.inverse();
    oldInvU = oldU.inverse();
    if (imp) {
        diff(oldPositions,qForceOld);
        qForceOld *= 2;
    } // end if
    newD = oldD;
    newU = oldU;
    while (true) {
        while (cycles < maxIterations) {
            /* run Monte Carlo cycles */
            for (unsigned int i = 0; i < oldPositions.rows(); ++i) {
                /* loop over number of particles(move only 1 particle) */
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

                // update (probability distribution function)
                if (i < halfSize) {
                    b->updateTrialWaveFunction(newD, newPositions.row(i),
                            alpha, i/2);
                } else {
                    b->updateTrialWaveFunction(newU, newPositions.row(i),
                            alpha, i/2);
                } // end ifelse

                if (imp) {
                    /* set new quantum force */
                    updateDiff(newPositions, qForceNew, i);
                    qForceNew.row(i) *= 2;
                } // end if

                // calculate Greens function ratio
                if (imp) {
                    transitionRatio = exp(0.125*step*(qForceOld.row(i).norm() -
                                qForceNew.row(i).norm()) +
                            0.25*step*((oldPositions(i,0)-newPositions(i,0)) *
                                (qForceNew(i,0)+qForceOld(i,0)) +
                                (oldPositions(i,1)-newPositions(i,1)) *
                                (qForceNew(i,1)+qForceOld(i,1))));
                } // end if

                if ((i<halfSize)) {
                    determinantRatioD = meth->determinantRatio(newD, oldInvD,
                            i/2);
                } else if ((i>=halfSize)) {
                    determinantRatioU = meth->determinantRatio(newU, oldInvU,
                            i/2);
                } // end ifelseif

                testRatio = determinantRatioD*determinantRatioD *
                    determinantRatioU*determinantRatioU * (!perturb ?  1 :
                            b->jastrowRatio(oldPositions, newPositions, beta,
                                i));
                if (imp) {
                    /* importance sampling */
                    testRatio *= transitionRatio;
                } //end if

                if (testRatio >= dist(mt)) {
                    /* update positions according to Metropolis test */
                    oldPositions.row(i) = newPositions.row(i);
                    if (i<halfSize) {
                        oldD.row(i/2) = newD.row(i/2);
                    } else {
                        oldU.row(i/2) = newU.row(i/2);
                    } // end if
                    if (imp) {
                        qForceOld.row(i) = qForceNew.row(i);
                    } // end if
                } else {
                    /* reset position */
                    newPositions.row(i) = oldPositions.row(i);
                } // end if

                // update inverse
                if (i < halfSize) {
                    meth->updateMatrixInverse(oldD, newD, oldInvD, newInvD,
                            determinantRatioD, i/2);
                } else {
                    meth->updateMatrixInverse(oldU, newU, oldInvU, newInvU,
                            determinantRatioU, i/2);
                } // end ifelse
            } // end fori

            // calculate local energy and local energy squared
            tmpEnergy = localEnergy2(oldPositions,perturb);
    //         tmpEnergy = localEnergyDiff(newD,newU,newPositions,perturb);
            energy += tmpEnergy;
            energySq += tmpEnergy*tmpEnergy;

            // calculate values for Hessen matrix
            tmpRsum = 0;
            tmphsum = 0;
            tmphsum3 = 0;
            tmpELalpRsum = 0;
            tmpELbethsum = 0;
            for (unsigned int k = 0; k < newPositions.rows(); ++k) {
                nkx = *(b->states[k][0]);
                nky = *(b->states[k][1]);
                tmpRsum -= newPositions.row(k).norm();
                Hxfactor = nkx*(nkx-1) *
                    H(newPositions(k,0),nkx-2)/H(newPositions(k,0),nkx);
                Hyfactor = nky*(nky-1) *
                    H(newPositions(k,1),nky-2)/H(newPositions(k,1),nky);
                ELalpRsum += 2*(nkx+nky+1) - Hxfactor - Hyfactor -
                    alpha*b->omega*newPositions.row(k).squaredNorm();
                for (unsigned int l = 0; l < newPositions.rows(); ++l) {
                    if (k != l) {
                        a = b->padejastrow(k,l);
                        rkl = (newPositions.row(k) -
                                newPositions.row(l)).norm();
                        denom = beta + 1 / rkl;
                        denomsq = denom*denom;
                        denomcu = denom*denom*denom;
                        tmphsum += a / denomsq;
                        tmphsum3 += a / denomcu;
                        tmpELalpRsum += 0.5*a*b->omega/denomcu *
                            (newPositions(k,0) *
                             (newPositions(k,0)-newPositions(l,0)) +
                             newPositions(k,1) *
                             (newPositions(k,1)-newPositions(l,1)));
                        tmpELbethsum += a/denomcu*(2 + 3*beta*beta*rkl +
                                2*rkl*(2*((Hxfactor-alpha*b->omega) *
                                        (newPositions(k,0)-newPositions(l,0)) +
                                        (Hyfactor-alpha*b->omega) *
                                        (newPositions(k,1)-newPositions(l,1)))
                                    - 1/denomsq));
                    } // end if
                } // end forj
            } // end fori
            Rsum += tmpRsum;
            hsum += tmphsum;
            Rsqsum += tmpRsum*tmpRsum;
            hsqsum += tmphsum*tmphsum;
            ELR += tmpEnergy*tmpRsum;
            ELh += tmpEnergy*tmphsum;
            ELRsq += tmpEnergy*tmpRsum*tmpRsum;
            ELhsq += tmpEnergy*tmphsum*tmphsum;
            ELRh += tmpEnergy*Rsum*hsum;
            ELalpRsum += tmpELalpRsum*tmpRsum;
            ELbethsum += tmpELbethsum*tmphsum;
            hsum3 += tmphsum3; 
            ELhsum3 += tmphsum3*tmpEnergy;

            cycles++;
        } // end while

        // calculate final energy estimation
        energy /= cycles;
        energySq /= cycles;
        Rsum /= cycles;
        hsum /= cycles;
        Rsqsum /= cycles;
        hsqsum /= cycles;
        hsum3 /= cycles;
        ELR /= cycles;
        ELh /= cycles;
        ELRsq /= cycles;
        ELhsq /= cycles;
        ELRh /= cycles;
        ELalpRsum /= cycles;
        ELbethsum /= cycles;
        ELhsum3 /= cycles;
        
        // set Hessen matrix
        HessenMatrix(0,0) = b->omega*b->omega*ELRsq - b->omega*ELalpRsum +
            0.5*b->omega*b->omega * (ELRsq + energy*Rsum*Rsum -
                    3*b->omega*energy*Rsum);
//         HessenMatrix(1,0) = 2*(b->omega*(ELRh - Rsum*ELh) -
//                 hsum*(energy*(2-b->omega*Rsum) + b->omega*ELR));
//         HessenMatrix(0,1) = HessenMatrix(1,0);
//         std::cout << HessenMatrix(0,1) << " " << HessenMatrix(1,0) << std::endl;
        HessenMatrix(0,1) = 0;
        HessenMatrix(1,0) = 0;
        HessenMatrix(1,1) = 4*ELhsq - 2*ELbethsum + 2*(ELhsq +
                2*energy*hsum*hsum) - 6*ELh*hsum + 4*(hsum3*energy + ELhsum3);

        // optimalize with CG
        rhs(0) = -energy;
        rhs(1) = -energy;
        newAlphaBeta = meth->conjugateGradient(HessenMatrix, rhs, newAlphaBeta);

        // conditional break
//         if (fabs(newAlphaBeta(0)-alpha)<=1e-10 &&
//                 fabs(newAlphaBeta(1)-beta)<=1e-10) {
//             /* break when variational parameters are steady */
//             break;
//         } // end if

        // set new variational parameters
        alpha = newAlphaBeta(0);
        beta = newAlphaBeta(1);
        std::cout << "alpha: " << alpha << ", beta: " << beta << std::endl;

        // Reset set variables used in Monte Carlo loop
        energy = 0;
        energySq = 0;
        Rsum = 0;
        hsum = 0;
        Rsqsum = 0;
        hsqsum = 0;
        hsum3 = 0;
        ELR = 0;
        ELh = 0;
        ELRsq = 0;
        ELhsq = 0;
        ELRh = 0;
        ELalpRsum = 0;
        ELbethsum = 0;
        ELhsum3 = 0;
        cycles = 0;
    } // end while true
} // end function calculate

//////////////////////////////////////////////////////////////////////////////
// Class for calculating variational Monte Carlo energies.                  //
//                                                                          //
// Main Functions:                                                          //
// See the individual functions for specific behavior.                      //
//////////////////////////////////////////////////////////////////////////////

#include "vmc.h" // header
#include <iomanip>
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

    aw = alpha*b->omega;
    awsqr = sqrt(aw);

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

void VMC::oneBodyFirstDerivativeRatio(const Eigen::MatrixXd &wave, const
        Eigen::MatrixXd &waveInv, Eigen::MatrixXd &der, const Eigen::MatrixXd
        &R, const unsigned int k, const unsigned int kIdx, const unsigned
        jstart) {
    /* Analytic first derivative ratio of one body part of wave function  for
     * particle k */
    double n;
    for (unsigned int d = 0; d < R.cols(); ++d) {
        n = *(b->states[jstart][d]);
        for (unsigned int i = 0; i < R.rows()/2; ++i) {
            der(k,d) += (2*awsqr*n*H(awsqr*R(k,d),n-1)/H(awsqr*R(k,d),n) -
                    aw*R(k,d)) * wave(kIdx,0)*waveInv(0,i);
        } // end forj
    } // end ford
} // end function oneBodyFirstDerivativeRatio

void VMC::oneBodySecondDerivativeRatio(const Eigen::MatrixXd &wave, const
        Eigen::MatrixXd &waveInv, Eigen::MatrixXd &der, const Eigen::MatrixXd
        &R, const unsigned int k, const unsigned int kIdx, const unsigned int
        jstart){
    /* Analytic second derivative of one body part of wave function for
     * particle k */
    double n, factor;
    for (unsigned int d = 0; d < R.cols(); ++d) {
        n = *(b->states[jstart][d]);
        factor = aw*(4*n*(n-1)*(n-1) * H(awsqr*R(k,d),n-2)/H(awsqr*R(k,d),n) +
                aw*R(k,d)*R(k,d) - 1 - n);
        for (unsigned int i = 0; i < R.rows()/2; ++i) {
            der(kIdx) += factor * wave(kIdx,0)*waveInv(0,i);
        } // end fori
    } // end ford
} // end function oneBodySecondDerivativeRatio

void VMC::jastrowFirstDerivativeRatio(Eigen::MatrixXd &der, const
        Eigen::MatrixXd &R, const int k) {
    /* Analytic first derivative ratio of Jastrow part of wave function for
     * particle k */
    double rjk;
    der.row(k).setZero();
    for (unsigned int j = 0; j < R.rows(); ++j) {
        if (j != k) {
            rjk = (R.row(k) - R.row(j)).norm();
            der.row(k) += b->padejastrow(k,j) * (R.row(k)-R.row(j)) /
                (rjk*pow(1+beta*rjk,2));
        } // end if
    } // end forj
} // end function jastrowFirstDerivativeRatio

double VMC::jastrowSecondDerivativeRatio(const Eigen::MatrixXd &R, const int k) {
    /* Analytic second derivative ratio of Jastrow part of wave function for
     * particle k */
    double ratio = 0;
    double rkj, rki, aki, akj, denomi, denomj;
    for (unsigned int i = 0; i < R.rows(); ++i) {
        if (i != k) {
            aki = b->padejastrow(k,i);
            rki = (R.row(k) - R.row(i)).norm();
            denomi = 1 + beta*rki;
            if (j != k) {
                for (unsigned int j = 0; j < R.rows(); ++j) {
                    akj = b->padejastrow(k,j);
                    rkj = (R.row(k) - R.row(j)).norm();
                    denomj = 1 + beta*rkj;
                    ratio += (R.row(k)-R.row(i)).dot(R.row(k)-R.row(j)) / (rki*rkj)
                        *aki*akj / (denomi*denomi*denomj*denomj);
                } // end forj
            } // end if
        } // end if
    } // end fori
    for (unsigned int j = 0; j < R.rows(); ++j) {
        if (j != k) {
            rkj = (R.row(k) - R.row(j)).norm();
            denomj = 1 + beta*rkj;
            ratio += 2*b->padejastrow(k,j) / (rkj*denomj*denomj*denomj);
        } // end if
    } // end forj
//     for (unsigned int j = 0; j < R.rows(); ++j) {
//         if (j != k) {
//             akj = b->padejastrow(k,j);
//             rkj = (R.row(k)-R.row(j)).norm();
//             denom = 1 + beta*rkj;
//             ratio += akj/(rkj*(denom*denom)) * (1 - 2*rkj*beta/denom);
//             for (unsigned int d = 0; d < R.cols(); ++d) {
//                 sumr(d) += akj * (R(k,d)-R(j,d)) / rkj*denom*denom;
//             } // end ford
//         } // end if
//     } // end forj
    return ratio;
} // end function jastrowSecondDerivativeRatio

double coulombFactor(const Eigen::MatrixXd &R) {
    // Calculate coulomb interaction */
    double C = 0;
    for (unsigned int i = 0; i < R.rows(); ++i) {
        for (unsigned int j = i+1; j < R.rows(); ++j) {
            C += 1 / (R.row(i)-R.row(j)).norm();
        } // end forj
    } // end fori
    return C;
} // end function coulombFactor

double VMC::localEnergy2(const Eigen::MatrixXd &lapD, const Eigen::MatrixXd
        &lapU, const Eigen::MatrixXd &derOB, const Eigen::MatrixXd &derJ, const
        Eigen::MatrixXd &R, bool coulomb) {
    /* calculate analytic expression of local energy */
    double E = 0;
    double halfSize = R.rows()/2;
    for (unsigned int k = 0; k < R.rows(); ++k) {
        E += 0.5*(b->omega*b->omega*R.row(k).squaredNorm() - (k<halfSize ?
                    lapD(k) : lapU(k-halfSize)) - (coulomb ?
                    jastrowSecondDerivativeRatio(R,k) -
                    2*derOB.row(k).dot(derJ.row(k)) : 0));
    } // end fork
    E += (coulomb ? coulombFactor(R) : 0);
    return E;
} // end function localEnergy2

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
    double testRatio, tmpEnergy, tmpR, tmpB2, tmpB3, tmpELalp, tmpELbet,
           transitionRatio, denom, denomsq, denomcu, a, rkl, Hxfactor,
           Hyfactor;
    unsigned int halfIdx;
    double A = 0;
    double RB2 = 0;
    double ELB3 = 0;
    double B3 = 0;
    double R = 0;
    double B2 = 0;
    double B2sq = 0;
    double Rsq = 0;
    double ELR = 0;
    double ELB2 = 0;
    double ELRsq = 0;
    double ELB2sq = 0;
    double ELRB2 = 0;
    double ELalpR = 0;
    double ELbetR = 0;
    double ELalpB2 = 0;
    double ELbetB2 = 0;
    double determinantRatioD = 1;
    double determinantRatioU = 1;
    unsigned int cycles = 0;
    unsigned int uIdx;
    int nkx, nky;
    double *determinantRatio;
    Eigen::MatrixXd derOB = Eigen::MatrixXd::Zero(oldPositions.rows(),
            oldPositions.cols());
    Eigen::MatrixXd derJ = Eigen::MatrixXd::Zero(oldPositions.rows(),
            oldPositions.cols());
    Eigen::MatrixXd lapU = Eigen::VectorXd::Zero(b->ECut);
    Eigen::MatrixXd lapD = Eigen::VectorXd::Zero(b->ECut);
    Eigen::MatrixXd *oldWave, *newWave, *oldInv, *newInv, *lap;
    Eigen::MatrixXd oldD = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd oldU = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd newD = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd newU = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd oldInvD = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd oldInvU = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd newInvD = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd newInvU = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    Eigen::MatrixXd HessenMatrix = Eigen::MatrixXd::Zero(2,2);
    Eigen::MatrixXd rhs = Eigen::VectorXd::Zero(2,1);
    Eigen::MatrixXd newAlphaBeta = Eigen::VectorXd::Zero(2,1);
    newAlphaBeta(0) = alpha;
    newAlphaBeta(1) = beta;
    b->setTrialWaveFunction(oldD, oldU, oldPositions, alpha);
    oldInvD = oldD.inverse();
    oldInvU = oldU.inverse();
    for (unsigned int k = 0; k < oldPositions.rows(); ++k) {
        /* set first derivatives */
        if (k<halfSize) {
            oneBodyFirstDerivativeRatio(oldD, oldInvD, derOB, oldPositions, k,
                    k, 0);
            oneBodySecondDerivativeRatio(oldD, oldInvD, lapD, oldPositions, k,
                    k, 0);
        } else {
            oneBodyFirstDerivativeRatio(oldU, oldInvU, derOB, oldPositions, k,
                    k-halfSize, 1);
            oneBodySecondDerivativeRatio(oldU, oldInvU, lapU, oldPositions, k,
                    k-halfSize, 1);
        } // end if
        jastrowFirstDerivativeRatio(derJ, oldPositions, k);
    } // end fork
    if (imp) {
        /* set quantum force */
        qForceOld = 2*(derOB + derJ);
    } // end if
    newD = oldD;
    newU = oldU;
    newInvD = oldInvD;
    newInvU = oldInvU;
    while (true) {
        for (cycles = 0; cycles < maxIterations; ++cycles) {
            /* run Monte Carlo cycles */
            for (unsigned int i = 0; i < oldPositions.rows(); ++i) {
                /* loop over number of particles(move only 1 particle) */
                // set references to matrices used (update only spin-up or
                // spin-down depending on which particle moved)
                if (i<halfSize) {
                    /* spin down for first N/2 particles */
                    oldWave = &oldD;
                    newWave = &newD;
                    oldInv = &oldInvD;
                    newInv = &newInvD;
                    determinantRatio = &determinantRatioD;
                    lap = &lapD;
                    halfIdx = i;
                    uIdx = 0;
                } else {
                    /* spin up for remaining N/2+1 to N particles */
                    oldWave = &oldU;
                    newWave = &newU;
                    oldInv = &oldInvU;
                    newInv = &newInvU;
                    determinantRatio = &determinantRatioU;
                    lap = &lapU;
                    halfIdx = i - halfSize;
                    uIdx = 1;
                } // end if
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

                // update Slater matrix
                b->updateTrialWaveFunction(*newWave, newPositions.row(i),
                        alpha, halfIdx, uIdx);

                // calculate determinant ratio
                *determinantRatio = meth->determinantRatio(*newWave, *oldInv,
                        halfIdx);

                // set newInv
                (*(newInv)).setZero();
                meth->updateMatrixInverse(*oldWave, *newWave, *oldInv, *newInv,
                        *determinantRatio, halfIdx);

                // update first derivatives
                if (imp) {
                    derOB.row(i).setZero();
                    derJ.row(i).setZero();
                    oneBodyFirstDerivativeRatio(*newWave, *newInv, derOB,
                            newPositions, i, halfIdx, uIdx);
                    jastrowFirstDerivativeRatio(derJ, newPositions, i);
                } // end if

                if (imp) {
                    /* set new quantum force */
                    qForceNew.row(i) = 2*(derOB.row(i) + derJ.row(i));
                } // end if

                // calculate transition function ratio for Metropolis test
                if (imp) {
                    transitionRatio = exp(0.125*step*(qForceOld.row(i).norm() -
                                qForceNew.row(i).norm()) +
                            0.25*step*((oldPositions(i,0)-newPositions(i,0)) *
                                (qForceNew(i,0)+qForceOld(i,0)) +
                                (oldPositions(i,1)-newPositions(i,1)) *
                                (qForceNew(i,1)+qForceOld(i,1))));
                } // end if

                // set Metropolis test
                testRatio = pow(*determinantRatio,2) * (!perturb ?  1 :
                        exp(2*b->jastrowRatio(oldPositions, newPositions, beta,
                                i)));

//                 if (i < halfSize) {
//                     std::cout << determinantRatioU << "U" << std::endl;
//                 } else {
//                     std::cout << determinantRatioD << "D" << std::endl;
//                 } // end if
                if (imp) {
                    /* importance sampling */
                    testRatio *= transitionRatio;
                } //end if

                if (testRatio >= dist(mt)) {
                    /* update state according to Metropolis test */
                    A++;
                    *oldInv = *newInv;
                    oldPositions.row(i) = newPositions.row(i);
                    oldWave->row(halfIdx) = newWave->row(halfIdx);
                    if (imp) {
                        qForceOld.row(i) = qForceNew.row(i);
                    } // end if
                } else {
                    /* reset(discard) state */
//                     *newInv = *oldInv;
//                     newPositions.row(i) = oldPositions.row(i);
//                     newWave->row(halfIdx) = oldWave->row(halfIdx);
//                     if (imp) {
//                         qForceNew.row(i) = qForceOld.row(i);
//                     } // end if
                } // end ifelse

                // update laplacian
                (*(lap))(halfIdx) = 0;
                oneBodySecondDerivativeRatio(*oldWave, *oldInv, *lap,
                        oldPositions, i, halfIdx, uIdx);
                *determinantRatio = meth->determinantRatio(*oldWave, *oldInv,
                        halfIdx);

                // update first derivatives
                if (imp) {
                    derOB.row(i).setZero();
                    derJ.row(i).setZero();
                    oneBodyFirstDerivativeRatio(*oldWave, *oldInv, derOB,
                            oldPositions, i, halfIdx, uIdx);
                    jastrowFirstDerivativeRatio(derJ, oldPositions, i);
                } // end if
            } // end fori

            // calculate local energy and local energy squared
            tmpEnergy = localEnergy2(lapD, lapU, derOB, derJ, oldPositions,
                    perturb); 
//             tmpEnergy = localEnergy2(oldPositions, perturb);
            energy += tmpEnergy;
            energySq += tmpEnergy*tmpEnergy;
                
            // calculate values for Hessen matrix
//             tmpR = 0;
//             tmpB2 = 0;
//             tmpB3 = 0;
//             tmpELalp = 0;
//             tmpELbet = 0;
//             for (unsigned int k = 0; k < newPositions.rows(); ++k) {
//                 nkx = *(b->states[k][0]);
//                 nky = *(b->states[k][1]);
//                 tmpR += newPositions.row(k).norm();
//                 Hxfactor = nkx*(nkx-1) *
//                     H(newPositions(k,0),nkx-2)/H(newPositions(k,0),nkx);
//                 Hyfactor = nky*(nky-1) *
//                     H(newPositions(k,1),nky-2)/H(newPositions(k,1),nky);
//                 tmpELalp += b->omega*(nkx+nky+1 + (Hxfactor + Hyfactor) -
//                         alpha*b->omega*newPositions.row(k).squaredNorm());
//                 for (unsigned int l = 0; l < newPositions.rows(); ++l) {
//                     if (k != l) {
//                         a = b->padejastrow(k,l);
//                         rkl = (newPositions.row(k) -
//                                 newPositions.row(l)).norm();
//                         denom = beta + 1 / rkl;
//                         denomsq = denom*denom;
//                         denomcu = denom*denom*denom;
//                         tmpB2 += a / denomsq;
//                         tmpB3 += a / denomcu;
//                         tmpELalp += a*b->omega/(rkl*denomsq) *
//                             (newPositions(k,0) *
//                              (newPositions(k,0)-newPositions(l,0)) +
//                              newPositions(k,1) *
//                              (newPositions(k,1)-newPositions(l,1)));
 //                         tmpELbet += a/denomcu * (rkl/denom *
 //                                 (2*beta+a*(1+beta)/denom - 1/rkl) +
 //                                 2*(((nkx+Hxfactor)/newPositions(k,0) -
 //                                         alpha*b->omega*newPositions(k,0)) *
 //                                     (newPositions(k,0)-newPositions(l,0)) +
 //                                     ((nky+Hyfactor)/newPositions(k,1) -
 //                                      alpha*b->omega*newPositions(k,1)) *
 //                                     (newPositions(k,1)-newPositions(l,1)))
 //                                 );
//                         tmpELbet -= a/denomcu * ((beta*rkl-2)/denom -
//                                 (newPositions(k,0)-newPositions(l,0) +
//                                  newPositions(k,1)-newPositions(l,1)) -
//                                 (((nkx+Hxfactor)/newPositions(k,0) -
//                                   alpha*b->omega*newPositions(k,0)) *
//                                  (newPositions(k,0)-newPositions(l,0)) +
//                                  ((nky+Hyfactor)/newPositions(k,1) -
//                                   alpha*b->omega*newPositions(k,1)) *
//                                  (newPositions(k,1)-newPositions(l,1))));
//                     } // end if
//                 } // end fork
//             } // end forl
//             R += tmpR;
//             B2 += tmpB2;
//             B3 += tmpB3; 
//             RB2 += tmpR*tmpB2;
//             Rsq += tmpR*tmpR;
//             B2sq += tmpB2*tmpB2;
//             ELB3 += tmpEnergy*B3;
//             ELR += tmpEnergy*tmpR;
//             ELB2 += tmpEnergy*tmpB2;
//             ELRB2 += tmpEnergy*R*B2;
//             ELalpR += tmpELalp*tmpR;
//             ELalpB2 += tmpELalp*tmpB2;
//             ELbetB2 += tmpELbet*tmpB2;
//             ELbetR += tmpELbet*tmpR;
//             ELRsq += tmpEnergy*tmpR*tmpR;
//             ELB2sq += tmpEnergy*tmpB2*tmpB2;
        } // end for cycles

        // calculate final expectation values
//         energy /= cycles * newPositions.rows();
//         energySq /= cycles * newPositions.rows();
        energy /= cycles;
        energySq /= cycles;

        std::cout << "Acceptance: " << A/(cycles*newPositions.rows()) << std::endl;
        break;
        R /= cycles;
        B2 /= cycles;
        RB2 /= cycles;
        Rsq /= cycles;
        B2sq /= cycles;
        B3 /= cycles;
        ELR /= cycles;
        ELB2 /= cycles;
        ELRsq /= cycles;
        ELB2sq /= cycles;
        ELRB2 /= cycles;
        ELalpR /= cycles;
        ELbetB2 /= cycles;
        ELB3 /= cycles;
        
        // set Hessen matrix
        HessenMatrix(0,0) = b->omega*(b->omega*(ELRsq - energy*Rsq + energy*R*R
                    - ELR*R) - ELalpR);
//         HessenMatrix(0,1) = 0;
//         HessenMatrix(1,0) = 0;
        HessenMatrix(0,1) = b->omega * (2*(ELRB2 - energy*RB2) +
                0.5*energy*R*B2 - ELR*B2 - ELB2*R - ELbetR);
        HessenMatrix(1,0) = b->omega * (2*(ELRB2 - energy*RB2) +
                0.5*energy*R*B2 - ELR*B2 - ELB2*R - ELalpB2);
        HessenMatrix(1,1) = 4*(ELB2sq - energy*B2sq) + energy*B2*B2 - ELB2*B2 -
            2*ELbetB2;

        std::cout << "Hessen: \n" << "  " << HessenMatrix(0,0) << " " << HessenMatrix(0,1) << "\n" << "  " << HessenMatrix(1,0) << " " << HessenMatrix(1,1) << std::endl;

        // optimalize with CG
        rhs(0) = b->omega*(energy*R - ELR);
        rhs(1) = 2*(energy*B2 - ELB2);
        newAlphaBeta += HessenMatrix.inverse()*rhs;

        // conditional break
//         if (fabs(newAlphaBeta(0)-alpha)<=1e-5 &&
//                 fabs(newAlphaBeta(1)-beta)<=1e-5) {
//             /* break when variational parameters are steady */
//             break;
//         } // end if

        // set new variational parameters
        std::cout << std::setprecision(10) << "alpha: " << alpha << ", beta: " << beta << ", Energy: " << energy << std::endl;
        alpha = newAlphaBeta(0);
        beta = newAlphaBeta(1);
       
        // stupid
        aw = alpha*b->omega;
        awsqr = sqrt(aw);

        // Reset variables used in Monte Carlo loop
        energy = 0;
        energySq = 0;
        R = 0;
        B2 = 0;
        Rsq = 0;
        B2sq = 0;
        B3 = 0;
        ELR = 0;
        ELB2 = 0;
        ELRsq = 0;
        ELB2sq = 0;
        ELRB2 = 0;
        ELalpR = 0;
        ELbetB2 = 0;
        ELB3 = 0;
    } // end while true
} // end function calculate

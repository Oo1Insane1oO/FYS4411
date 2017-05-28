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
#include <sstream>
#include <iterator>
#include <time.h>
#include <string.h>
#include <cstdio>

VMC::VMC(Basis *B, double alp, double bet, unsigned int d, double s, unsigned
        int max) {
    alpha = alp;
    beta = bet;
    b = B;
    dim = d;
    step = s;
    maxIterations = max;

    aw = alpha*b->omega;
    awsqr = sqrt(aw);

    imp = false;
    coulomb = false;
    jastrow = false;

    meth = new Methods(); 
} // end constructor

VMC::~VMC() {
    delete meth;
} // end deconstructor

void VMC::setImportanceSampling(bool a) {
    /* switch importance sampling on */
    imp = a;
} //send function setImportanceSampling

void VMC::setCoulombInteraction(bool a) {
    /* switch Coloumb interaction on */
    coulomb = a;
} // end function setCoulombInteraction

void VMC::setJastrow(bool a) {
    /* switch Jastrow factor on */
    jastrow = a;
} // end function setCoulombInteraction

void VMC::setAllOn() {
    /* switch importance sampling, Coulomb interaction and Jastrow factor on */
    setImportanceSampling(true);
    setCoulombInteraction(true);
    setJastrow(true);
} // end function setAllOn

void VMC::oneBodyFirstDerivativeRatio(const Eigen::MatrixXd &wave, const
        Eigen::MatrixXd &waveInv, Eigen::MatrixXd &der, const Eigen::MatrixXd
        &R, const unsigned int k, const unsigned int kIdx, const unsigned
        jstart) {
    /* Analytic first derivative ratio of one body part of wave function  for
     * particle k */
    int n;
    for (unsigned int j = 0; j < R.rows(); j+=2) {
        for (unsigned int d = 0; d < R.cols(); ++d) {
            n = *(b->states[j+jstart][d]);
            der(k,d) += (2*awsqr*n*H(awsqr*R(k,d),n-1)/H(awsqr*R(k,d),n) -
                    aw*R(k,d)) * wave(kIdx,j/2) * waveInv(j/2,kIdx);
        } // end ford
    } // end forj
} // end function oneBodyFirstDerivativeRatio

void VMC::oneBodySecondDerivativeRatio(const Eigen::MatrixXd &wave, const
        Eigen::MatrixXd &waveInv, Eigen::MatrixXd &der, const Eigen::MatrixXd
        &R, const unsigned int k, const unsigned int kIdx, const unsigned int
        jstart){
    /* Analytic second derivative of one body part of wave function for
     * particle k */
    for (unsigned int j = 0; j < R.rows(); j+=2) {
        for (unsigned int d = 0; d < R.cols(); ++d) {
            der(kIdx) += aw * (aw*R(k,d)*R(k,d) - 1 -
                    2**(b->states[j+jstart][d])) * wave(kIdx,j/2)
                * waveInv(j/2,kIdx);
        } // end ford
    } // end forj
} // end function oneBodySecondDerivativeRatio

void VMC::jastrowFirstDerivativeRatio(Eigen::MatrixXd &der, const
        Eigen::MatrixXd &R, const unsigned int k) {
    /* Analytic first derivative ratio of Jastrow part of wave function for
     * particle k */
    double rjk;
    for (unsigned int j = 0; j < R.rows(); ++j) {
        if (j != k) {
            rjk = (R.row(k) - R.row(j)).norm();
            der.row(k) += b->padejastrow(k,j) * (R.row(k)-R.row(j)) /
                (rjk*pow(1+beta*rjk,2));
        } // end if
    } // end forj
} // end function jastrowFirstDerivativeRatio

double VMC::jastrowSecondDerivativeRatio(const Eigen::MatrixXd &R, const
        unsigned int k) {
    /* Analytic second derivative ratio of Jastrow part of wave function for
     * particle k */
    double ratio = 0;
    double rkj, denom;
    for (unsigned int j = 0; j < R.rows(); ++j) {
        if (j != k) {
            rkj = (R.row(k) - R.row(j)).norm();
            denom = 1 + beta*rkj;
            ratio += b->padejastrow(k,j)/(rkj*denom*denom) * (dim - 1 -
                    2*beta*rkj/denom);
        } // end if
    } // end forj
    return ratio + derJ.row(k).squaredNorm();
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
        Eigen::MatrixXd &R) {
    /* calculate analytic expression of local energy */
    double E = 0;
    double halfSize = R.rows()/2;
    for (unsigned int k = 0; k < R.rows(); ++k) {
        E += 0.5*(b->omega*b->omega*R.row(k).squaredNorm() - (k<halfSize ?
                    lapD(k) : lapU(k-halfSize)) - (jastrow ?
                    jastrowSecondDerivativeRatio(R,k) +
                    2*derOB.row(k).dot(derJ.row(k)) : 0));
    } // end fork
    return (coulomb ? (E + coulombFactor(R)) : E);
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

void VMC::setSeed(long int s) {
    /* set seed */
    seed = s;
} // end function setSeed

long int VMC::getSeed() {
    /* get seed */
    return seed;
} // end function setSeed

double VMC::Afunc(const Eigen::MatrixXd &wave, const Eigen::MatrixXd &waveInv,
        const Eigen::MatrixXd &R, const unsigned int iStart) {
    /* first derivative of wave function with respect to alpha */
    double A = 0;
    double n;
    for (unsigned int d = 0; d < R.cols(); ++d) {
        for (unsigned int i = 0; i < 2*R.rows(); i+=2) {
                n = *(b->states[i+iStart][d]);
            for (unsigned int j = 0; j < R.rows(); ++j) {
                A += n/(alpha) * (sqrt(alpha) + 2*R(j,d)*(n-1)*sqrt(b->omega) *
                        H(awsqr*R(j,d),n-2)/H(awsqr*R(j,d),n) -
                        b->omega*R.row(j).squaredNorm()) * wave(j,i/2) *
                    waveInv(i/2,j);
             } // end forj
        } // end fori
    } // end ford
    return A;
} // end function Afunc

double VMC::Bfunc(const Eigen::MatrixXd &R) {
    /* first derivative of wave function with respect to beta */
    double B = 0;
    for (unsigned int i = 0; i < R.rows(); ++i) {
        for (unsigned int j = 0; j < R.rows(); ++j) {
            if (i != j)  {
                B -= b->padejastrow(i,j) / pow((beta +
                            1/(R.row(i)-R.row(j)).norm()),2);
            } // end if
        } // end forj
    } // end fori
    return B;
} // end function Bfunc

void VMC::initializeCalculationVariables() {
    /* initialize Eigen matrices/vectors used in function calculate */
    oldPositions = Eigen::MatrixXd::Zero(2*b->ECut, dim);
    newPositions = Eigen::MatrixXd::Zero(2*b->ECut, dim);

    // force vector in steepest descent
    steepb = Eigen::VectorXd::Zero(2);
    prevSteepb = Eigen::VectorXd::Zero(2);

    // first derivatives
    derOB = Eigen::MatrixXd::Zero(oldPositions.rows(),
            oldPositions.cols());
    derJ = Eigen::MatrixXd::Zero(oldPositions.rows(), oldPositions.cols());
   
    // spin down matrices
    lapD = Eigen::VectorXd::Zero(b->ECut);
    oldD = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    newD = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    oldInvD = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    newInvD = Eigen::MatrixXd::Zero(b->ECut, b->ECut);

    // spin up matrices
    lapU = Eigen::VectorXd::Zero(b->ECut);
    oldU = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    newU = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    oldInvU = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    newInvU = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
   
    // variational parameters
    newAlphaBeta = Eigen::VectorXd::Zero(2);
    oldAlphaBeta = Eigen::VectorXd::Zero(2);
        
    if (imp) {
        qForceOld = Eigen::MatrixXd::Zero(oldPositions.rows(),
                oldPositions.cols());
        qForceNew = Eigen::MatrixXd::Zero(oldPositions.rows(),
                oldPositions.cols());
    } // end if
} // end function initializeCalculationVariables

void VMC::calculate(const char *destination) {
    /* function for running Monte Carlo integration */

    // initialize Mersenne Twister random number generator and uniform
    // distribution engine
    unsigned int halfIdx, uIdx;
    double testRatio, tmpEnergy, acceptance, tmpA, tmpB, A, ELA, B, ELB;

    double determinantRatioD = 1;
    double determinantRatioU = 1;
    unsigned int cycles = 0;
    double *determinantRatio;
        
    Eigen::MatrixXd *oldWave, *newWave, *oldInv, *newInv, *lap;
   
    std::ofstream openFile;

    // initialize random number generator
    std::istringstream stringBuffer("0 1 2 3 4 5 6 7 8 9 10");
    std::istream_iterator<int> start(stringBuffer), end;
    std::seed_seq seedSequence(start, end);
//     std::mt19937_64 mt(time(NULL));
    std::mt19937_64 mt(seedSequence);
    std::uniform_real_distribution<double> dist(0,1);
    std::normal_distribution<double> normDist(0,1);

    // set sizes
    initializeCalculationVariables();
    unsigned int halfSize = oldPositions.rows()/2;
    double steepStep = 0.01;

    unsigned int runCount = 1;

    char tmpf[80];

    while (runCount <= 500) {
        // reinitialize positions
        for (unsigned int i = 0; i < oldPositions.rows(); ++i) {
            for (unsigned int j = 0; j < oldPositions.cols(); ++j) {
                if (imp) {
                    oldPositions(i,j) = normDist(mt) * sqrt(step);
                } else {
                    oldPositions(i,j) = step * (dist(mt)-0.5);
                } // end ifelse
            } // end forj
        } // end fori

        // set variational parameter vector
        oldAlphaBeta(0) = alpha;
        oldAlphaBeta(1) = beta;
        newAlphaBeta = oldAlphaBeta;

        // initialize Slater matrix and its inverse
        b->setTrialWaveFunction(oldD, oldU, oldPositions, alpha);
        oldInvD = oldD.inverse();
        oldInvU = oldU.inverse();

        for (unsigned int k = 0; k < oldPositions.rows(); ++k) {
            /* set first derivatives */
            if (k<halfSize) {
                if (jastrow || imp) {
                    oneBodyFirstDerivativeRatio(oldD, oldInvD, derOB, oldPositions,
                            k, k, 0);
                } // end if
                oneBodySecondDerivativeRatio(oldD, oldInvD, lapD, oldPositions, k,
                        k, 0);
            } else {
                if (jastrow || imp) {
                    oneBodyFirstDerivativeRatio(oldU, oldInvU, derOB, oldPositions,
                            k, k-halfSize, 1);
                } // end if
                oneBodySecondDerivativeRatio(oldU, oldInvU, lapU, oldPositions, k,
                        k-halfSize, 1);
            } // end if
            if (jastrow) {
                jastrowFirstDerivativeRatio(derJ, oldPositions, k);
            } // end if
        } // end fork

        if (imp) {
            /* set quantum force */
            qForceOld = 2*(derOB + derJ);
        } // end if
        
        // reset values used in Monte Carlo cycle
        energy = 0;
        energySq = 0;
        A = 0;
        ELA = 0;
        B = 0;
        ELB = 0;
        acceptance = 0;
    
        if (destination) {
            sprintf(tmpf, "%s_%d.txt", destination, runCount);
            openFile.open(tmpf);
        } // end ifelseif

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

                // update Slater, ratio and inverse
                b->updateTrialWaveFunction(*newWave, newPositions, alpha, i,
                        halfIdx, uIdx);
                *determinantRatio = meth->determinantRatio(*newWave, *oldInv,
                        halfIdx);
                (*(newInv)).setZero();
                meth->updateMatrixInverse(*oldWave, *newWave, *oldInv, *newInv,
                        *determinantRatio, halfIdx);

                // update first derivatives
                if (imp) {
                    derOB.row(i).setZero();
                    oneBodyFirstDerivativeRatio(*newWave, *newInv, derOB,
                            newPositions, i, halfIdx, uIdx);
                } // end if

                if (jastrow) {
                    derJ.row(i).setZero();
                    jastrowFirstDerivativeRatio(derJ, newPositions, i);
                } // end if

                if (imp) {
                    /* set new quantum force */
                    qForceNew.row(i) = 2*(derOB.row(i) + derJ.row(i));
                } // end if

                // set Metropolis test
                testRatio = *determinantRatio * *determinantRatio * (jastrow ?
                        exp(2*b->jastrowRatio(oldPositions, newPositions, beta,
                                i)) : 1);
                if (imp) {
                    /* importance sampling, calculate transition function ratio
                     * for Metropolis test */
                    testRatio *= exp(0.125*step*(qForceOld.row(i).norm() -
                                qForceNew.row(i).norm()) +
                            0.25*step*((oldPositions(i,0)-newPositions(i,0)) *
                                (qForceNew(i,0)+qForceOld(i,0)) +
                                (oldPositions(i,1)-newPositions(i,1)) *
                                (qForceNew(i,1)+qForceOld(i,1))));
                } //end if

                if (testRatio >= dist(mt)) {
                    /* update state according to Metropolis test */
                    acceptance++;
                    *oldInv = *newInv;
                    oldPositions.row(i) = newPositions.row(i);
                    oldWave->row(halfIdx) = newWave->row(halfIdx);
                    if (imp) {
                        qForceOld.row(i) = qForceNew.row(i);
                    } // end if
                } // end if

                // update Laplacian and determinant ratio
                (*(lap))(halfIdx) = 0;
                oneBodySecondDerivativeRatio(*oldWave, *oldInv, *lap,
                        oldPositions, i, halfIdx, uIdx);
                *determinantRatio = meth->determinantRatio(*oldWave, *oldInv,
                        halfIdx);

                // update first derivatives
                if (imp) {
                    derOB.row(i).setZero();
                    oneBodyFirstDerivativeRatio(*oldWave, *oldInv, derOB,
                            oldPositions, i, halfIdx, uIdx);
                } // end if
                if (jastrow) {
                    derJ.row(i).setZero();
                    jastrowFirstDerivativeRatio(derJ, oldPositions, i);
                } // end if
            } // end fori

            // calculate local energy and local energy squared
            tmpEnergy = localEnergy2(lapD, lapU, derOB, derJ, oldPositions);
            openFile << tmpEnergy << " " << tmpEnergy*tmpEnergy << "\n";
            energy += tmpEnergy;
            energySq += tmpEnergy*tmpEnergy;

            // split spin up/down and calculate expected value(local) of first
            // derivative of wave function with respect to alpha
            tmpA = Afunc(oldD, oldInvD, oldPositions.block(0, 0, halfSize,
                        oldPositions.cols()), 0) + Afunc(oldU, oldInvU,
                    oldPositions.block(halfSize, 0, halfSize,
                        oldPositions.cols()),1);
            A += tmpA;
            ELA += tmpEnergy*tmpA;

            // No need for splitting when finding first derivative with respect
            // to beta(only Jastrow factor gives constribution)
            tmpB = Bfunc(oldPositions);
            B += tmpB;
            ELB += tmpEnergy*tmpB;
        } // end for cycles

        // calculate final expectation values
        energy /= cycles;
        energySq /= cycles;
        A /= cycles;
        ELA /= cycles;
        ELB /= cycles;
        B /= cycles;

//         openFile << " " << "\n";
//         openFile << "Summed total: " << energy << " " << energySq << "\n";
//         openFile << "Acceptance: " << acceptance/(cycles*newPositions.rows()) << "\n";
//         openFile << "Alpha: " << alpha << "\n";
//         openFile << "Beta: " << beta << "\n";
//         openFile.close();
        openFile << " " << "\n";
        openFile << energy << " " << energySq << "\n";
        openFile << acceptance/(cycles*newPositions.rows()) << "\n";
        openFile << alpha << "\n";
        openFile << beta << "\n";
        openFile.close();

        std::cout << "Acceptance: " << acceptance/(cycles*newPositions.rows()) << std::endl;

        // optimalize with steepest descent method
        steepb(0) = ELA - energy*A;
        steepb(1) = 2*(ELB - energy*B);
        newAlphaBeta = oldAlphaBeta - steepStep*steepb;

        std::cout << "alpha: " << alpha << " beta: " << beta << " Energy: " <<
            energy << std::endl;

        // update stepsize in steepest descent according to two-step size
        // gradient
        steepStep = (newAlphaBeta.row(0) -
                oldAlphaBeta.row(0)).transpose().dot(steepb.row(0) -
                prevSteepb.row(0)) / (steepb - prevSteepb).squaredNorm();

        // update variational parameters
        alpha = newAlphaBeta(0);
        beta = newAlphaBeta(1);
        prevSteepb = steepb;
        oldAlphaBeta = newAlphaBeta;
        aw = alpha*b->omega;
        awsqr = sqrt(aw);
      
        runCount++;
    } // end while true
} // end function calculate

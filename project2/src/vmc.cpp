//////////////////////////////////////////////////////////////////////////////
// Class for calculating variational Monte Carlo energies.                  //
//                                                                          //
// Main Functions:                                                          //
// See the individual functions for specific behavior.                      //
//////////////////////////////////////////////////////////////////////////////

#include "vmc.h" // header
#include "hermite.h" // template function for hermite polynomials
#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include <sstream>
#include <iterator>
#include <time.h>
#include <string.h>
#include <cstdio>
#include <cmath>
#include <cstdlib>

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
    
//     std::istringstream stringBuffer("0 1 2 3 4 5 6 7 8 9 10");
//     std::istream_iterator<int> start(stringBuffer), end;
//     std::seed_seq seedSequence(start, end);
//     mt.seed(seedSequence);
    std::mt19937_64 mt(time(NULL));
    dist = std::uniform_real_distribution<double>(0,1);
    normDist = std::normal_distribution<double>(0,1);
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

void VMC::oneBodyFirstDerivativeRatio(Eigen::MatrixXd &buf, const
        Eigen::MatrixXd &R, const unsigned int k, const unsigned int j) {
    /* Analytic first derivative ratio of one body part of wave function  for
     * particle k */
    int n;
    for (unsigned int d = 0; d < R.cols(); ++d) {
        n = *(b->states[j][d]);
        buf(0,d) = (2*awsqr*n*H(awsqr*R(k,d),n-1)/H(awsqr*R(k,d),n) -
                aw*R(k,d));
    } // end ford
} // end function oneBodyFirstDerivativeRatio

double VMC::oneBodySecondDerivativeRatio(const Eigen::MatrixXd &R, const
        unsigned int k, const unsigned int j) {
    double sum = 0;
    for (unsigned int d = 0; d < dim; ++d) {
        sum += aw * (aw*R(k,d)*R(k,d) - 1 - 2**(b->states[j][d]));
    } // end ford
    return sum;
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
    jbuf.setZero();
    jastrowFirstDerivativeRatio(jbuf,R,k);
    return ratio + jbuf.row(k).squaredNorm();
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

double VMC::calculateLocalEnergy(const Eigen::MatrixXd &waveD, const
        Eigen::MatrixXd &waveU, const Eigen::MatrixXd &waveInvD, const
        Eigen::MatrixXd &waveInvU, const Eigen::MatrixXd &R, const
        Eigen::MatrixXd &derOB, const Eigen::MatrixXd &derJ) {
    /* Analytic expression for local energy */
    double E = 0;
    unsigned int half = R.rows()/2;
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Zero(1,dim);
    for (unsigned int k = 0; k < R.rows(); ++k) {
        /* Potential part */
        E += 0.5*b->omega*b->omega*R.row(k).squaredNorm();
        for (unsigned int j = 0; j < R.rows(); j+=2) {
            /* One-body Laplacian part*/
            if (k < half) {
                /* spin down */
                E -= 0.5 * oneBodySecondDerivativeRatio(R,k,j) * waveD(k,j/2) *
                    waveInvD(j/2,k);
            } else {
                /* spin up */
                E -= 0.5 * oneBodySecondDerivativeRatio(R,k,j+1) *
                    waveU(k-half,j/2) * waveInvU(j/2,k-half);
            } // end ifelse
        } // end forj
    } // end fork
    if (jastrow) {
        /* set Jastrow part */
        for (unsigned int k = 0; k < R.rows(); ++k) {
            E -= 0.5*jastrowSecondDerivativeRatio(R,k);
            tmp.setZero();
            for (unsigned int j = 0; j < R.rows(); j+=2) {
                if (k<half) {
                    oneBodyFirstDerivativeRatio(buf,R,k,j);
                    tmp += buf * waveD(k,j/2) * waveInvD(j/2,k);
                } else {
                    oneBodyFirstDerivativeRatio(buf,R,k,j+1);
                    tmp += buf * waveU(k-half,j/2) * waveInvU(j/2,k-half);
                } // end if
            } // end forj
            E -= tmp.row(0).dot(jbuf.row(k));
        } // end fork
    } // end if
    return (coulomb ? E + coulombFactor(R) : E);
} // end if

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

double VMC::Afunc(const Eigen::MatrixXd &wave, const Eigen::MatrixXd &waveInv,
        const Eigen::MatrixXd &R, const unsigned int iStart) {
    /* first derivative of wave function with respect to alpha */
    double A = 0;
    double n;
    for (unsigned int i = 0; i < 2*R.rows(); i+=2) {
        for (unsigned int j = 0; j < R.rows(); ++j) {
            A -= 0.5*b->omega*R.row(j).squaredNorm();
            for (unsigned int d = 0; d < R.cols(); ++d) {
                n = *(b->states[i+iStart][d]);
                A += 1/(2*alpha*sqrt(alpha)) * n*(1+R(j,d)*(n-1) *
                        sqrt(b->omega/alpha) *
                        H(awsqr*R(j,d),n-2)/H(awsqr*R(j,d),n));
            } // end ford
        } // end forj
    } // end fori
    return A;
} // end function Afunc

double VMC::Bfunc(const Eigen::MatrixXd &R) {
    /* first derivative of wave function with respect to beta */
    double B = 0;
    for (unsigned int i = 0; i < R.rows(); ++i) {
        for (unsigned int j = 0; j < R.rows(); ++j) {
            if (i != j)  {
                B -= b->padejastrow(i,j) / pow(beta +
                        1/(R.row(i)-R.row(j)).norm(),2);
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
    derOB = Eigen::MatrixXd::Zero(oldPositions.rows(), oldPositions.cols());
    derJ = Eigen::MatrixXd::Zero(oldPositions.rows(), oldPositions.cols());
   
    // spin down matrices
    oldD = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    newD = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    oldInvD = Eigen::MatrixXd::Zero(b->ECut, b->ECut);
    newInvD = Eigen::MatrixXd::Zero(b->ECut, b->ECut);

    // spin up matrices
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
    
    buf = Eigen::MatrixXd(1,dim);
    jbuf = Eigen::MatrixXd(oldPositions.rows(),dim);
} // end function initializeCalculationVariables

void VMC::setFirstDerivatives(const Eigen::MatrixXd &wave, const
        Eigen::MatrixXd &waveInv, const Eigen::MatrixXd &R, const unsigned int
        k, const unsigned int kIdx, const unsigned int jstart) {
    /* Calculate and set first derivatives for particle k */
    derOB.row(k).setZero();
    for (unsigned int j = 0; j < R.rows(); j+=2) {
        oneBodyFirstDerivativeRatio(buf, R, k, j+jstart);
        derOB.row(k) += buf * wave(kIdx,j/2) * waveInv(j/2,kIdx);
    } // end forj
    if (jastrow) {
        derJ.row(k).setZero();
        jastrowFirstDerivativeRatio(derJ,R,k);
    } // end if
} // end function setFirstDerivatives 
        
void VMC::initializePositions(Eigen::MatrixXd &R) {
    /* Set initial positions */
    for (unsigned int i = 0; i < R.rows(); ++i) {
        for (unsigned int j = 0; j < R.cols(); ++j) {
            if (imp) {
                R(i,j) = normDist(mt) * sqrt(step);
            } else {
                R(i,j) = step * (dist(mt)-0.5);
            } // end ifelse
        } // end forj
    } // end fori
} // end function initializePositions

void VMC::calculate(const char *destination) {
    /* function for running Monte Carlo integration */
    unsigned int halfIdx, uIdx, randomDim;
    double testRatio, tmpEnergy, acceptance, tmpA, tmpB, A, ELA, B, ELB;

    double determinantRatioD = 1;
    double determinantRatioU = 1;
    unsigned int cycles = 0;
    double *determinantRatio;
        
    Eigen::MatrixXd *oldWave, *newWave, *oldInv, *newInv;
   
    std::ofstream openFile;
    
    // set sizes
    initializeCalculationVariables();
    unsigned int halfSize = oldPositions.rows()/2;
    double steepStep = 0.01;

    // count for filenames and buffer for filename string
    unsigned int runCount = 1;
    char tmpf[80];

    while (runCount <= 1000) {
        // reinitialize positions
        initializePositions(oldPositions);
        newPositions = oldPositions;

        // set variational parameter vector
        oldAlphaBeta(0) = alpha;
        oldAlphaBeta(1) = beta;
        newAlphaBeta = oldAlphaBeta;

        // initialize Slater matrix and its inverse
        b->setTrialWaveFunction(oldD, oldU, oldPositions, alpha);
        oldInvD = oldD.inverse();
        oldInvU = oldU.inverse();
        newD = oldD;
        newU = oldU;

        if (imp || jastrow) {
            for (unsigned int k = 0; k < oldPositions.rows(); ++k) {
                if (k<halfSize) {
                    setFirstDerivatives(oldD, oldD, oldPositions, k, k, 0);
                } else {
                    setFirstDerivatives(oldU, oldU, oldPositions, k,
                            k-halfSize, 1);
                } // end if
            } // end fork
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

        unsigned int i;
        for (cycles = 0; cycles < maxIterations; ++cycles) {
            /* run Monte Carlo cycles */
            /* loop over number of particles(move only 1 particle) */
            i = rand() % oldPositions.rows();

            // set references to matrices used (update only spin-up or
            // spin-down depending on which particle moved)
            if (i<halfSize) {
                /* spin down for first N/2 particles */
                oldWave = &oldD;
                newWave = &newD;
                oldInv = &oldInvD;
                newInv = &newInvD;
                determinantRatio = &determinantRatioD;
                halfIdx = i;
                uIdx = 0;
            } else {
                /* spin up for remaining N/2+1 to N particles */
                oldWave = &oldU;
                newWave = &newU;
                oldInv = &oldInvU;
                newInv = &newInvU;
                determinantRatio = &determinantRatioU;
                halfIdx = i - halfSize;
                uIdx = 1;
            } // end if
            
            // propose new position and move only in one dimension
            randomDim = (std::fabs(dist(mt)) < 0.5 ? 0 : 1);
            if (imp) {
                newPositions(i,randomDim) = oldPositions(i,randomDim) +
                    0.5*qForceOld(i,randomDim)*step + normDist(mt)*sqrt(step);
            } else {
                newPositions(i,randomDim) = oldPositions(i,randomDim) +
                    step*(dist(mt)-0.5);
            } // end ifelse

            // update Slater matrix, determinant ratio and Slater inverse
            b->updateTrialWaveFunction(*newWave, newPositions, alpha, i,
                    halfIdx, uIdx);
            *determinantRatio = meth->determinantRatio(*newWave, *oldInv,
                    halfIdx);
            newInv->setZero();
            meth->updateMatrixInverse(*oldWave, *newWave, *oldInv, *newInv,
                    *determinantRatio, halfIdx);

            // update first derivatives (ratios)
            if (imp || jastrow) {
                setFirstDerivatives(*newWave, *newInv, newPositions, i,
                        halfIdx, uIdx);
            } // end if

            if (imp) {
                /* set new quantum force */
                qForceNew.row(i) = 2*(derOB.row(i) + derJ.row(i));
            } // end if

            // set Metropolis test
            testRatio = *determinantRatio * *determinantRatio * (jastrow ?
                    exp(2*b->jastrowRatio(oldPositions, newPositions, beta,
                            i)) : 1.0);
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

            if (testRatio >= dist(mt) || testRatio > 1) {
                /* update state according to Metropolis test */
                acceptance++;
                *oldInv = *newInv;
                oldPositions.row(i) = newPositions.row(i);
                oldWave->row(halfIdx) = newWave->row(halfIdx);
                if (imp) {
                    qForceOld.row(i) = qForceNew.row(i);
                } // end if
            } else {
                /* reset(discard) proposed state */
                *newInv = *oldInv;
                newPositions.row(i) = oldPositions.row(i);
                newWave->row(halfIdx) = oldWave->row(halfIdx);
                if (imp) {
                    qForceOld.row(i) = qForceNew.row(i);
                } // end if
            } // end if

            // update determinant ratio
            *determinantRatio = meth->determinantRatio(*oldWave, *oldInv,
                    halfIdx);

            // update first derivatives
            if (imp || jastrow) {
                setFirstDerivatives(*oldWave, *oldInv, oldPositions, i,
                        halfIdx, uIdx);
            } // end if
       
            // Accumulate local energy and local energy squared
            tmpEnergy = calculateLocalEnergy(oldD, oldU, oldInvD, oldInvU,
                    oldPositions, derOB, derJ);
            energy += tmpEnergy;
            energySq += tmpEnergy*tmpEnergy;

            // write to file
            if (destination) {
                openFile << tmpEnergy << " " << tmpEnergy*tmpEnergy << "\n";
            } // end if

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

        if (destination) {
            openFile << " " << "\n";
            openFile << energy << " " << energySq << "\n";
            openFile << acceptance/(cycles*newPositions.rows()) << "\n";
            openFile << alpha << "\n";
            openFile << beta << "\n";
            openFile.close();
        } // end if

        std::cout << "Acceptance: " << acceptance/cycles << std::endl;

        // optimalize with steepest descent method
        steepb(0) = 2*(ELA - energy*A);
        steepb(1) = 2*(ELB - energy*B);
        newAlphaBeta = oldAlphaBeta - steepStep*steepb;

        std::cout << std::setprecision(10) << "alpha: " << alpha << " beta: "
            << beta << " Energy: " << energy << std::endl;

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

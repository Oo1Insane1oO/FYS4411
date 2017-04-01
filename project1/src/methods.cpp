//////////////////////////////////////////////////////////////////////////////
// Class containing functions for different methods.                        //
//                                                                          //
// Functions:                                                               //
//      - Kronk: function for Kronecker delta.                              //
//      - factorial: recursive function for positive integer factorial.     //
//      - findRootsHermite: find roots of a polynomial.                     //
//      - gaussHermiteQuadrature: solve integral of integrand of form       //
//        exp(-x^2)H_n(x)^2.                                                //
//      - generateRandom: Generate a pair of number based on Box-Muller.    //
// See the individual functions for specific behavior.                      //
//////////////////////////////////////////////////////////////////////////////

#include "methods.h"
#include "Coulomb_Functions.hpp"
#include <math.h>
#include <iostream>
#include <random>

Methods::Methods() {
} // end constructor

int Methods::kronk(int *i, int *j) {
    /* function for Kroncecker delta function delta = 1 if i=j, 0 otherwise */
    return ((i==j) ? 1 : 0);
} // end function kronk

int Methods::kronk(int &i, int &j) {
    /* function for Kroncecker delta function delta = 1 if i=j, 0 otherwise */
    return ((i==j) ? 1 : 0);
} // end function kronk

int Methods::factorial(int x) {
    /* recursively find factorial */
    return (x==1 ? x : (x==0 ? 1 : x*factorial(x-1)));
} // end function factorial

double Methods::hermiteNormal(int n) {
    /* normalization constant */
    return 1./(pow(2,n) * factorial(n) * sqrt(M_PI));
} // end function hermiteNormal

double Methods::hermite(double x, int n) {
    /* calculate hermite polynomial of degree n */
    if(n<0) {
        /* negative indices dont exist */
        return 0;
    } else if(n==0) {
        /* first value, H0 */
        return 1;
    } else {
        /* recursive relation */
        return 2*x*hermite(x,n-1) - 2*(n-1)*hermite(x,n-2);
    } // endifelseifelse
} // end function hermite

void Methods::copyVector(std::vector<double>& v, std::vector<double> w) {
    /* copy vector w onto vector v */
    for (unsigned int i = 0; i < v.size(); ++i) {
        v[i] = w[i];
    } // end fori
} // end function copyVector

void Methods::findRootsHermite(int level, std::vector<double> &x) {
    /* find roots of hermite polynomial of order level with Aberth-Erlich
     * method */

    // variables for later use
    int size = x.size();
    double tmpSum = 0;
    double dx = 0.001;
    double dxHalf = dx/2;

    // copy original vector
    std::vector<double> prevP(level);
    copyVector(prevP,x);

    // run algorithm
    bool diffTest = true;
    while(diffTest) {
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if(i!=j) {
                    tmpSum += 1/(x[i]-x[j]);
                } // end if j not i
            } // end forj
            prevP[i] += 1/(tmpSum - (hermite(x[i]+dxHalf,level) -
                        hermite(x[i]-dxHalf,level)) /
                    (dx*hermite(x[i],level)));
            tmpSum = 0;
        } // end fori
        for (unsigned int k = 0; k < prevP.size(); ++k) {
            if(fabs(x[k]-prevP[k]) <= 1e-14) {
                /* check that all roots are within convergence */
                prevP[k] = (fabs(prevP[k])<=1e-14) ? 0 : prevP[k];
                diffTest = false;
            } else {
                /* break and rerun if convergence is not reached */
                diffTest = true;
                break;
            } // end ifelse
        } // end fork
        copyVector(x,prevP);
    } // end while diffTest
} // end function findRoots

void Methods::setPoints(int size, std::vector<double> &x) {
    /* set points for integral.
     * Points are the roots of the hermite polynomial of order size */
    x.resize(size,0);

    // set initial guess for roots (random uniform distribution)
    double sizeInv = 1./size;
    std::vector<double> pair(2);
    for (int i = 0; i < size; ++i) {
        generateRandom(pair);
        x[i] = pair[0]*size + pair[1]*sizeInv;
    } // end fori

    // find roots
    findRootsHermite(size,x);
} // end function set points

void Methods::setWeight(int size, std::vector<double> &weight, std::vector<double> &x) {
    /* set the weights for Gauss-Hermite quadrature */
    // set points, allocate weights vector, set a for later use
    int s = size<2 ? 2 : size;
    weight.resize(s,0);
    setPoints(s,x);
    double a = pow(2,s-1)*factorial(s)*sqrt(M_PI)/(s*s);
    double h;
    for (int i = 0; i < s; ++i) {
        h = hermite(x[i],s-1);
        weight[i] = a/(h*h);
    } // end fori
} // end function setWeights

int Methods::checkn(int varn) {
    return ((varn==0 || varn==1) ? 2 : varn);
} // end function checkn

double Methods::gaussHermiteQuadrature(int *n,
        std::function<double(double,double,double,double)> V) {
    /* Find integral based on given points n 
     * Assume V and H_(n_i)(x_i) commutes */
    double hi, hj, hk, hl;
    std::vector<std::vector<double>> p(4);
    std::vector<std::vector<double>> w(4);
    double A = 1.;
    for (int i = 0; i < 4; ++i) {
        /* set points, weights and normalization constant */
        setWeight(n[i]+1,w[i],p[i]);
        A *= hermiteNormal(n[i]);
    } // end fori

    double sum = 0;
    for (unsigned int i = 0; i < p[0].size(); ++i) {
        hi = hermite(p[0][i],n[0]);
        for (unsigned int j = 0; j < p[1].size(); ++j) {
            hj = hermite(p[1][j],n[1]);
            for (unsigned int k = 0; k < p[2].size(); ++k) {
                hk = hermite(p[2][k],n[2]);
                for (unsigned int l = 0; l < p[3].size(); ++l) {
                    hl = hermite(p[3][l],n[3]);
                    sum += w[0][i]*w[1][j]*w[2][k]*w[3][l] * A *
                        hi*hi*hj*hj*hk*hk*hl*hl *
                        V(p[0][i],p[1][j],p[2][k],p[3][l]);
                } // end forl
            } // end fork
        } // end forj
    } // end fori
    return sum;
} // end function gaussHermiteQuadrature

void Methods::generateRandom(std::vector<double> &pair) {
    /* generate pair of random numbers based on Box-Muller transform */
    double u,v,s,S;
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> randomReal(-1,1);
    while (true) {
        u = randomReal(generator);
        v = randomReal(generator);
        s = u*u + v*v;
        if(s>0 && s<1) {
            break;
        } // endif s
    } // end while true
    S = sqrt(-2*log(s)/s);
    pair[0] = S*u;
    pair[1] = S*v;
} // end function generateRandom

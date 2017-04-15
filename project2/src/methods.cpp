//////////////////////////////////////////////////////////////////////////////
// Class containing functions for different methods.                        //
//                                                                          //
// Functions:                                                               //
//      - factorial: recursive function for positive integer factorial.     //
// See the individual functions for specific behavior.                      //
//////////////////////////////////////////////////////////////////////////////

#include "methods.h"
#include <math.h>
#include <iostream>

Methods::Methods() {
} // end constructor

int Methods::factorial(int x) {
    /* recursively find factorial */
    return (x==1 ? x : (x==0 ? 1 : x*factorial(x-1)));
} // end function factorial

double Methods::min(double var1, double var2) {
    return (var1<var2 ? var1 : var2);
} // end function min

int Methods::min(int var1, int var2) {
    return (var1<var2 ? var1 : var2);
} // end template function min

double Methods::max(double var1, double var2) {
    return (var1>var2 ? var1 : var2);
} // end function max

int Methods::max(int var1, int var2) {
    return (var1>var2 ? var1 : var2);
} // end template function max

void Methods::updateMatrixInverse(Eigen::MatrixXd Mold, Eigen::MatrixXd Mnew,
        Eigen::MatrixXd MoldInv, Eigen::MatrixXd &MnewInv, unsigned int i) {
    /* update inverse of matrix when only column i has changed */
    double R = determinantRatio(Mnew, MoldInv, i);
    unsigned int N = MnewInv.rows();
    for (unsigned int k = 0; k < MnewInv.rows(); ++k) {
        for (unsigned int j = 0; j < MnewInv.rows(); ++j) {
            for (unsigned int l = 0; l < MnewInv.rows(); ++l) {
                if (j==i) {
                    MnewInv(k,j) += Mold(i,l)*MoldInv(l,j);
                } else {
                    MnewInv(k,j) -= Mnew(i,l)*MoldInv(l,j);
                } // end ifelse
            } // end forl
            MnewInv(k,j) *= MoldInv(k,i)/R;
            MnewInv(k,j) = ((j!=i) ? MnewInv(k,j)+MoldInv(k,j) : MnewInv(k,j));
        } // end forj
    } // end fork
} // end function updateMatrixInverse

double Methods::determinantRatio(Eigen::MatrixXd newElement, Eigen::MatrixXd
        oldInverse, unsigned int i) {
    /* Calculate determinant ratio of Slater determinants */
    double R = 0;
    for (unsigned int j = 0; j < oldInverse.rows(); ++j) {
        R += newElement(i,j) * oldInverse(j,i);
    } // end fori
    return R;
} // end function determinantRatio

double Methods::variance(double p,double psq) {
    /* calculate variance given <p> and <p^2>, that is expectation value and
     * expectation value squared */
    return psq - p*p;
} // end function variance

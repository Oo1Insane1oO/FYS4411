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

void Methods::updateMatrixInverse(Eigen::MatrixXd Mold, Eigen::MatrixXd Mnew,
        Eigen::MatrixXd MoldInv, Eigen::MatrixXd &MnewInv, unsigned int i) {
    /* update inverse of matrix when only column i has changed */
    double R = determinantRatio(Mnew, Mold, i);
    unsigned int N = MnewInv.size();
    for (unsigned int k = 0; k < N; ++k) {
        for (unsigned int j = 0; j < N; ++j) {
            for (unsigned int l = 0; l < N; ++l) {
                if (j==i) {
                    MnewInv(k,j) += Mold(i,l)*MoldInv(i,j);
                } else {
                    MnewInv(k,j) -= Mnew(i,l)*MoldInv(i,j);
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
        R += newElement(j,i) * oldInverse(j,i);
    } // end fori
    return R;
} // end function determinantRatio

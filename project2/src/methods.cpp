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

Eigen::MatrixXd updateMatrixInverse(Eigen::MatrixXd rcurr, 
        Eigen::MatrixXd rnew, Eigen::MatrixXd &M) {
    /* update inverse of matrix when only one column has changed */
} // end function updateMatrixInverse

double Basis::determinantRatio(Eigen::MatrixXd currPos, 
        Eigen::MatrixXd newPos) {
    /* Calculate determinant ratio of Slater determinants */
} // end function determinantRatio

///////////////////////////////////////////////////////////////////////////////
// Class containing functions for different methods.                         //
//                                                                           //
// Functions:                                                                //
//      - factorial: recursive function for positive integer factorial.      //
//      - min: minimum value of given arguments.                             //
//      - max: maximum value of given arguments.                             //
//      - updateMatrixInverse: Update inverse of matrix with Shermin-Morrison//
//      formula.                                                             //
//      - determinantRatio: Calculate determinant ratio when only one row    //
//      differs.                                                             //
//      - variance: calculate and return variance of given arguments.        //
//      - conjugateGradient: Find x for equation Ax=b using the conjugate    //
//      gradient method.                                                     //
// See the individual functions for specific behavior.                       //
///////////////////////////////////////////////////////////////////////////////

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
    /* find and return minimum value */
    return (var1<var2 ? var1 : var2);
} // end function min

int Methods::min(int var1, int var2) {
    /* find and return minimum value */
    return (var1<var2 ? var1 : var2);
} // end template function min

double Methods::max(double var1, double var2) {
    /* find and return maximum value */
    return (var1>var2 ? var1 : var2);
} // end function max

int Methods::max(int var1, int var2) {
    /* find and return maximum value */
    return (var1>var2 ? var1 : var2);
} // end template function max

void Methods::updateMatrixInverse(const Eigen::MatrixXd &Mold, const
        Eigen::MatrixXd &Mnew, const Eigen::MatrixXd &MoldInv, Eigen::MatrixXd
        &MnewInv, const double &ratio, unsigned int i) {
    /* update inverse of matrix when only row i has changed */
    MnewInv.setZero();
    for (unsigned int k = 0; k < MnewInv.rows(); ++k) {
        for (unsigned int j = 0; j < MnewInv.rows(); ++j) {
            for (unsigned int l = 0; l < MnewInv.rows(); ++l) {
                if (j==i) {
                    MnewInv(k,j) += Mold(i,l)*MoldInv(l,j);
                } else {
                    MnewInv(k,j) -= Mnew(i,l)*MoldInv(l,j);
                } // end ifelse
            } // end forl
            MnewInv(k,j) *= MoldInv(k,i)/ratio;
            MnewInv(k,j) = ((j!=i) ? MnewInv(k,j)+MoldInv(k,j) : MnewInv(k,j));
        } // end forj
    } // end fork
} // end function updateMatrixInverse

double Methods::determinantRatio(const Eigen::MatrixXd &newElement, const
        Eigen::MatrixXd &oldInverse, unsigned int i) {
    /* Calculate determinant ratio of determinants */
    double R = 0;
    for (unsigned int j = 0; j < oldInverse.rows(); ++j) {
        R += newElement(i,j) * oldInverse(j,i);
    } // end fori
    return R;
} // end function determinantRatio

double Methods::variance(double p,double psq, int N) {
    /* calculate variance given <p> and <p^2>, that is expectation value and
     * expectation value squared */
    return 1./N * (psq - p*p);
} // end function variance

Eigen::MatrixXd Methods::conjugateGradient(const Eigen::MatrixXd &A, const
        Eigen::MatrixXd &rhs, const Eigen::MatrixXd &x0) {
    /* solve a linear system, Ax=b with Conjugate Gradient method */
    Eigen::MatrixXd res = rhs - A*x0;
    Eigen::MatrixXd p = res;
    Eigen::MatrixXd xold = x0;
    Eigen::MatrixXd xnew;
    Eigen::MatrixXd rnew = Eigen::MatrixXd::Zero(x0.rows(),x0.cols());
    double C = res.squaredNorm();
    double rInner;
    rnew << 1, 1;
    unsigned int iter = 0;
    while (rnew.norm() > 1e-14) {
        rInner = res.squaredNorm();
        C = rInner / (p.transpose()*A*p).sum();
        xnew = xold + C*p;
        rnew = res - C*A*p;
        p = rnew + rnew.squaredNorm()/rInner * p;
        res = rnew;
        xold = xnew;
        iter++;
    } // end while
    return xnew;
} // end function conjugateGradient

#include "hermite1.h"
#include <cmath>

double HermitePolynomials::evaluate(int n, double x, double omega, double alpha) {

    if(n == 0) {
        return 1.0;
    } else if(n == 1) {
        return 2.0*sqrt(omega*alpha)*x;
    } else if(n == 2) {
        return 4.0*omega*alpha*x*x - 2.0;
    } else if(n == 3) {
        double sqrt_omega_alpha = sqrt(omega*alpha);
        return 8.0*sqrt_omega_alpha*omega*alpha*x*x*x - 12.0*sqrt_omega_alpha*x;
    }
}

double HermitePolynomials::evaluateDerivative(int n, double x, double omega, double alpha) {
    if(n == 0) {
        return 0.0;
    } else if(n == 1) {
        return 2.0*sqrt(omega*alpha);
    } else if(n == 2) {
        return 8.0*x*omega*alpha;
    } else if(n == 3) {
        double sqrt_omega_alpha = sqrt(omega*alpha);
        return 24.0*sqrt_omega_alpha*omega*alpha*x*x - 12.0*sqrt_omega_alpha;
    }
}

double HermitePolynomials::evaluateDoubleDerivative(int n, double x, double omega, double alpha) {
    if(n == 0) {
        return 0.0;
    } else if(n == 1) {
        return 0.0;
    } else if(n == 2) {
        return 8.0*omega*alpha;
    } else if(n == 3) {
        return 48.0*sqrt(omega*alpha)*omega*alpha*x;
    } else if(n == 4) {
        return 192.0*x*x-96.0;
    }
}

double HermitePolynomials::evaluateAlphaDerivative(int n, double x, double omega, double alpha) {
    if(n == 0) {
        return 0.0;
    } else if(n == 1) {
        return x*sqrt(omega/alpha);
    } else if(n==2) {
        return 4*omega*x*x;
    } else if(n==3) {
        return 12*sqrt(omega*alpha)*omega*x*x*x - 6*x*sqrt(omega/alpha);
    }
}

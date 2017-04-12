//////////////////////////////////////////////////////////////////////////////
// Header file for class Basis, see .cpp file for more information          //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#ifndef BASIS_H
#define BASIS_H

#include "methods.h"
#include <Eigen/Dense>

class Basis {
    private:
        int ECut, s;
        std::vector<int> n, ms, E, M, m;
        double omega, alpha, beta, a;
        Eigen::MatrixXd phiU, phiD;

        Methods *meth;

        void pushState(std::vector<int*>&, int, int, int);
        double jastrow(double, double);

    public:
        Basis(double, int);
        virtual ~Basis();

        std::vector<std::vector<int*>> states;
        void printStates();

        void setBasisMatrix(Eigen::MatrixXd, double); 
        double harmonicOscillatorWaveFunction(double, double, int, int);
        double trialWaveFunction(Eigen::MatrixXd, double, double, double);
};

#endif /* BASIS_H */

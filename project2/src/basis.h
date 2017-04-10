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
        int ECut;
        double omega;

        Methods *meth;

    public:
        Basis(double, int);
        virtual ~Basis();

        int s;
        std::vector<int> n, ms, E, M, m;
        std::vector<std::vector<int*>> states;
        void printStates();

        double harmonicOscillatorWaveFunction(double, double, int, int);
        double trialWaveFunction(Eigen::MatrixXd, double, double, double);
};

#endif /* BASIS_H */

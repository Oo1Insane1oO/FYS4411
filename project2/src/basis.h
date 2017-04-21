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
        int s;
        std::vector<int> n, ms, E, M, m;
        Eigen::MatrixXd phi;

        Methods *meth;

        void pushState(std::vector<int*>&, int, int, int);
        void setSpinMatrix();

    public:
        Basis(double, int);
        virtual ~Basis();

        int ECut;
        double omega;
        std::vector<std::vector<int*>> states;

        double harmonicOscillatorWaveFunction(double, double, double, int,
                int);
        Eigen::MatrixXd trialWaveFunction(const Eigen::MatrixXd&, double);
        double jastrow(const Eigen::MatrixXd&, double);
        void printStates();
};

#endif /* BASIS_H */

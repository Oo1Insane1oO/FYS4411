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
        double jastrow(double, double, double, double);
        void setSpinMatrix();

    public:
        Basis(double, int);
        virtual ~Basis();

        int ECut;
        double omega;
        std::vector<std::vector<int*>> states;

        void setBasisMatrix(Eigen::MatrixXd, double); 
        double harmonicOscillatorWaveFunction(double, double, double, int,
                int);
        Eigen::MatrixXd trialWaveFunction(Eigen::MatrixXd, double, double);
        void printStates();
};

#endif /* BASIS_H */

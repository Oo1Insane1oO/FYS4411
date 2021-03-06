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

        Methods *meth;

        void pushState(std::vector<int*>&, int, int, int, int);
        void findPossiblenxny(unsigned int, std::vector<std::vector<int>>& p);

    public:
        Basis(double, unsigned int);
        virtual ~Basis();

        unsigned int ECut;
        double omega;
        std::vector<std::vector<int*>> states;

        double harmonicOscillatorWaveFunction(double, double, double, int,
                int);
        void setTrialWaveFunction(Eigen::MatrixXd&, Eigen::MatrixXd&, const
                Eigen::MatrixXd&, const double);
        void updateTrialWaveFunction(Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const double, const unsigned int, const unsigned int, const
                unsigned int);
        double jastrow(const Eigen::MatrixXd&, double);
        void updateJastrow(double&, const Eigen::MatrixXd&, const
                Eigen::MatrixXd&, double, unsigned int);
        double jastrowRatio(const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                double, unsigned int);
        double padejastrow(const unsigned int&, const unsigned int&);
        void printStates();

        std::vector<int> getMagicNumbers();
};

#endif /* BASIS_H */

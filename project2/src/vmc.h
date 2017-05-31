//////////////////////////////////////////////////////////////////////////////
// Header file for class vmc, see .cpp file for more information            //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#ifndef VMC_H
#define VMC_H

#include "basis.h"
#include <string>

class VMC {
    private:
        Basis *b;
        Methods *meth;

        long int seed = -1;

        double aw, awsqr;

        bool imp, coulomb, jastrow;

        double Afunc(const Eigen::MatrixXd&);
        double Afunc(const Eigen::MatrixXd&, const Eigen::MatrixXd&, const
                Eigen::MatrixXd&, const unsigned int jstart);
        double Bfunc(const Eigen::MatrixXd&);

        void setFirstDerivatives(const Eigen::MatrixXd&, const
                Eigen::MatrixXd&, const Eigen::MatrixXd&, const unsigned int,
                const unsigned int, const unsigned int);
        void initializePositions(Eigen::MatrixXd&);
        void initializeCalculationVariables();

        Eigen::MatrixXd oldPositions, newPositions, qForceOld, qForceNew,
            steepb, prevSteepb, derOB, derJ, oldD, oldU, newD, newU, oldInvD,
            oldInvU, newInvD, newInvU;

        Eigen::MatrixXd buf, jbuf;

        std::mt19937_64 mt;
        std::uniform_real_distribution<double> dist;
        std::normal_distribution<double> normDist;
    public:
        VMC (Basis*, double, double, unsigned int, double, unsigned int);
        virtual ~VMC ();

        unsigned int dim, maxIterations;
        double alpha, beta, energy, energySq, step;
        Eigen::MatrixXd newAlphaBeta, oldAlphaBeta;

        void diff(const Eigen::MatrixXd&, Eigen::MatrixXd&);
        void updateDiff(const Eigen::MatrixXd&, Eigen::MatrixXd&, unsigned
                int);
        double diff2(Eigen::MatrixXd&, Eigen::MatrixXd&, const
                Eigen::MatrixXd&, double);

        double localEnergy2(const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXd&, const Eigen::MatrixXd&, const
                Eigen::MatrixXd&);
        double localEnergy2(const Eigen::MatrixXd &R, bool=true);
        double localEnergyDiff(Eigen::MatrixXd&, Eigen::MatrixXd&, const
                Eigen::MatrixXd&, bool=true);

        void oneBodyFirstDerivativeRatio(Eigen::MatrixXd&, const
                Eigen::MatrixXd&, const unsigned int, const unsigned int);
        double oneBodySecondDerivativeRatio(const Eigen::MatrixXd&, const unsigned
                int, const unsigned int);
        double calculateLocalEnergy(const Eigen::MatrixXd&, const
                Eigen::MatrixXd&, const Eigen::MatrixXd&, const
                Eigen::MatrixXd&, const Eigen::MatrixXd&, const
                Eigen::MatrixXd&, const Eigen::MatrixXd&); 

        void jastrowFirstDerivativeRatio(Eigen::MatrixXd&, const
                Eigen::MatrixXd&, const unsigned int);
        double jastrowSecondDerivativeRatio(const Eigen::MatrixXd&, const unsigned int);

        void setImportanceSampling(bool);
        void setCoulombInteraction(bool);
        void setJastrow(bool);
        void setAllOn();

        void calculate(const unsigned int, const char* = NULL);
};

#endif /* VMC_H */

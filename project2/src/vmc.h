//////////////////////////////////////////////////////////////////////////////
// Header file for class vmc, see .cpp file for more information            //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#ifndef VMC_H
#define VMC_H

#include "basis.h"

class VMC {
    private:
        Basis *b;
        Methods *meth;

        unsigned long int seed = 1286754;

        double aw, awsqr;

    public:
        VMC (Basis*, double, double, unsigned int, double, unsigned int, bool);
        virtual ~VMC ();

        bool imp;
        unsigned int dim, maxIterations;
        double alpha, beta, energy, energySq, step;

        unsigned long int getSeed();
        void setSeed(unsigned long int);

        void diff(const Eigen::MatrixXd&, Eigen::MatrixXd&);
        void updateDiff(const Eigen::MatrixXd&, Eigen::MatrixXd&, unsigned
                int);
        double diff2(Eigen::MatrixXd&, Eigen::MatrixXd&, const
                Eigen::MatrixXd&, double);

        double localEnergy2(const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXd&, const Eigen::MatrixXd&, const
                Eigen::MatrixXd&, bool=true);
        double localEnergyDiff(Eigen::MatrixXd&, Eigen::MatrixXd&, const
                Eigen::MatrixXd&, bool=true);

        void oneBodyFirstDerivativeRatio(const Eigen::MatrixXd&, const
                Eigen::MatrixXd&, Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const unsigned int, const unsigned int, const unsigned int);
        void oneBodySecondDerivativeRatio(const Eigen::MatrixXd&, const
                Eigen::MatrixXd&, Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const unsigned int, const unsigned int);

        void jastrowFirstDerivativeRatio(Eigen::MatrixXd&, const
                Eigen::MatrixXd&, const int);
        double jastrowSecondDerivativeRatio(const Eigen::MatrixXd&, const int);

        void calculate(bool=true);
};

#endif /* VMC_H */

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

    public:
        VMC (Basis*, double, double, unsigned int, double, unsigned int);
        virtual ~VMC ();

        unsigned int dim, maxIterations;
        double alpha, beta, a, energy, energySq, step;

        double diff2(Eigen::MatrixXd, double);
        double localEnergy2(Eigen::MatrixXd, bool=true);
        double localEnergyDiff(Eigen::MatrixXd, bool=true);
        void calculate(unsigned long int=86754);
};

#endif /* VMC_H */

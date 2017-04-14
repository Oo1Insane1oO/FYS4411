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
        VMC (Basis*, double, double, unsigned int);
        virtual ~VMC ();

        unsigned int dim;
        double alpha, beta, a, energy, energySq;

        double diff2(Eigen::MatrixXd, double);
        double localEnergy2(Eigen::MatrixXd, bool=true);
        double localEnergyDiff(Eigen::MatrixXd, bool=true);
        void calculate(double, int, unsigned long int=86754);
};

#endif /* VMC_H */

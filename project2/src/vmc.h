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
        VMC (Basis*, double, double, double);
        virtual ~VMC ();

        double alpha, beta, a, energy, step;

        Eigen::MatrixXd R;

        double localEnergy2(Eigen::MatrixXd, Eigen::MatrixXd, bool=true);
        void initialize(long int=85754);
        void calculate();
        double metropolisTest(double, double);
};

#endif /* VMC_H */

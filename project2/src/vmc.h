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
        VMC (Basis*);
        virtual ~VMC ();

        double alpha, beta, a, energy;

        double localEnergy2(Eigen::MatrixXd, Eigen::MatrixXd, bool=true);
        double metropolisTest(double, double);
};

#endif /* VMC_H */

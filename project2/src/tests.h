#ifndef TESTS_H
#define TESTS_H

#include "vmc.h"

class Tests {
    private:
        Basis *b;
        VMC *v;
        Methods *m;

        Eigen::MatrixXd oldM, newM;
        int rowi;
        double eps;

        bool test_energies();
        bool test_energy();
        bool test_determinantratio();
        bool test_updateinverse();
        bool test_wavefunction2();
        bool test_updateWaveFunction();
        bool test_padejastrow();
        bool test_conjugateGradient();
        bool test_negativeHermite();
        bool test_hermite();
    
    public:
        Tests (Basis*, VMC*, int=3);
        virtual ~Tests ();

        void run_tests(int);
};

#endif /* TESTS_H */

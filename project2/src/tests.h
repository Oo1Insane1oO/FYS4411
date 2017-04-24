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
        bool test_2particle();
        bool test_determinantratio();
        bool test_updateinverse();
        bool test_wavefunction2();
    
    public:
        Tests (Basis*, VMC*, int=3);
        virtual ~Tests ();

        void run_tests(int);
};

#endif /* TESTS_H */

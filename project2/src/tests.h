#ifndef TESTS_H
#define TESTS_H

#include "vmc.h"

class Tests {
    private:
        Basis *b;
        VMC *v;
        Methods *m;

        bool test_energies();
        bool test_2particle();
    
    public:
        Tests (Basis*, VMC*);
        virtual ~Tests ();

        void run_tests(int);
};

#endif /* TESTS_H */

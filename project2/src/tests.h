#ifndef TESTS_H
#define TESTS_H

#include "basis.h"

class Tests {
    private:
        Basis *b;

        bool test_energies();
    
    public:
        Tests (Basis*);
        virtual ~Tests ();

        void run_tests(int);
};

#endif /* TESTS_H */

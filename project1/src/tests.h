//////////////////////////////////////////////////////////////////////////////
// Header file for class Tests                                              //
//////////////////////////////////////////////////////////////////////////////

#ifndef TESTS_H
#define TESTS_H

#include "basis.h"

class Tests {
    private:
        int A;
        Basis *b;
        Methods *meth;
        bool test_hermite();
        bool test_convert();
        bool test_energies();
        bool test_integrals();
        bool test_hartreefock();

    public:
        Tests (double, int, int);
        virtual ~Tests ();

        void run_tests(int);
};

#endif /* TESTS_H */

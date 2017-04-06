//////////////////////////////////////////////////////////////////////////////
// Class containing functions for different methods.                        //
//                                                                          //
// Functions:                                                               //
//      - factorial: recursive function for positive integer factorial.     //
// See the individual functions for specific behavior.                      //
//////////////////////////////////////////////////////////////////////////////

#include "methods.h"
#include <math.h>
#include <iostream>

Methods::Methods() {
} // end constructor

int Methods::factorial(int x) {
    /* recursively find factorial */
    return (x==1 ? x : (x==0 ? 1 : x*factorial(x-1)));
} // end function factorial

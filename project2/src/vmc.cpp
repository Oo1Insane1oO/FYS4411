//////////////////////////////////////////////////////////////////////////////
// Class for calculating variational Monte Carlo energies.                  //
//                                                                          //
// Main Functions:                                                          //
// See the individual functions for specific behavior.                      //
//////////////////////////////////////////////////////////////////////////////

#include "vmc.h" // header
#include <iostream>
#include <fstream>
#include <math.h>

VMC::VMC(Basis *B) {
    b = B;
} // end constructor

VMC::~VMC() {
    delete b;
} // end deconstructor

//////////////////////////////////////////////////////////////////////////////
// Class for assembling a harmonic oscillator basis for calculating         //
// variatonal monte Carlo energies.                                         //
//                                                                          //
// Main Functions:                                                          //
// See the individual functions for specific behavior.                      //
//                                                                          //
// Constructor sets up the states in the basis as a vector states. Each     //
// elements in vector states is a vector with 6 elements:                   //
//      nx: Principal quantum number(used for energy)                       //
//      ny: Principal quantum number(used for energy)                       //
//      s : Spin quantum number (always a whole number regardless of        //
//          fermion/boson                                                   //
//      ms: Spin projection quantum number (also a whole number)            //
//      E : Energy level given as nx+ny+1                                   //
//      M : Magic number(number of particles present at the given energy    // 
//          level of state                                                  //
//////////////////////////////////////////////////////////////////////////////

#include "basis.h" // header
#include <iostream>
#include <fstream>
#include <math.h>
#include <Eigen/Eigenvalues>

Basis::Basis(double w, int cut) {
    /* initialize states */
   
    // set omega
    omega = w;

    // set cutoff, spin and spin projections
    ECut = cut;
    s = 1;
    ms.resize(2);
    ms[0] = -1;
    ms[1] = 1;

    // set possible values for nx and ny
    n.resize(cut);
    E.resize(cut);
    for (unsigned int i = 0; i < n.size(); ++i) {
        n[i] = i;
        E[i] = i+1;
    } // end fori

    // allocate state arrays as {nx,ny,s,ms}
    M.resize(cut,0);
    std::vector<int*> s1 = std::vector<int*>(6,0);
    std::vector<int*> s2 = std::vector<int*>(6,0);
    for (int i = 0; i < ECut; ++i) {
        /* loop over values for nx */
        for (int j = 0; j <= i; ++j) {
            /* set states not yet pushed */
            for (int k = 0; k <= 1; ++k) {
                /* set both spin projections */
                if((i==j && k==1) || (i+j >= ECut)) {
                    /* dont count doubly and make sure cutoff is reached */
                    break;
                } // endif

                // set nx, ny and spin
                s1[0] = &(n[i]);
                s2[0] = &(n[j]);
                s1[1] = &(n[j]);
                s2[1] = &(n[i]);
                s1[2] = &s;
                s2[2] = &s;
              
                // set spin projection
                if(k==0) {
                    s1[3] = &(ms[k]);
                    s2[3] = &(ms[k+1]);
                } else if(k==1) {
                    s1[3] = &(ms[k]);
                    s2[3] = &(ms[k-1]);
                } // end ifelseif

                // set energy
                s1[4] = &(E[i+j]);
                s2[4] = &(E[i+j]);

                // increment and magic number
                M[i+j] += 1;
                s1[5] = &(M[i+j]);
                s2[5] = &(M[i+j]);
               
                // push state to states array
                states.push_back(s1);
                states.push_back(s2);
            } // end fork
        } // end forj
    } // end fori

    // sum magic numbers
    for (int i = M.size()-1; i >= 0; --i) {
        for (int j = 0; j <= i; ++j) {
            M[i] += M[j];
        } // end forj
    } // end fori
} // end constructor

Basis::~Basis() {
} // end deconstructor

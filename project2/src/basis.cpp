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
#include "hermite.h" // hermite polynomials
#include <iostream>
#include <fstream>

Basis::Basis(double w, int cut) {
    /* initialize states */
    
    // set methods object
    meth = new Methods();
   
    // set omega
    omega = w;

    // set cutoff, spin and spin projections
    ECut = cut;
    s = 1;
    ms.resize(2);
    ms[0] = -1;
    ms[1] = 1;

    // set possible values for nx and ny and energy E
    n.resize(cut);
    E.resize(cut);
    for (unsigned int i = 0; i < n.size(); ++i) {
        n[i] = i;
        E[i] = i+1;
    } // end fori

    // allocate state arrays as {nx,ny,s,ms,E,M}
    M.resize(cut,0);
    std::vector<int*> s1 = std::vector<int*>(6,0);
    for (int i = 0; i < ECut; ++i) {
        /* loop over values for nx */
        for (int j = 0; j <= i; ++j) {
            /* set states not yet pushed */
            
            if (i+j>=ECut) {
                /* end when cutoff is reached */
                break;
            } // end if

            // increment magic number and set values and push to states
            M[i+j]++;
            pushState(s1, i, j, 0);
            pushState(s1, i, j, 1);

            // dont set states doubly
            if (i!=j) {
                M[i+j]++;
                pushState(s1, j, i, 0);
                pushState(s1, j, i, 1);
            } // end if
        } // end forj
    } // end fori

    // sum magic numbers
    for (int i = M.size()-1; i >= 0; --i) {
        for (int j = 0; j <= i; ++j) {
            M[i] += M[j];
        } // end forj
    } // end fori
} // end constructor

void Basis::pushState(std::vector<int*> &state, int i, int j, int ud) {
    /* set (nx,ny,s,ms,E,M) in s1 and push s1 to states */
    state[0] = &(n[i]);
    state[1] = &(n[j]);
    state[2] = &s;
    state[3] = &(ms[ud]);
    state[4] = &(E[i+j]);
    state[5] = &(M[i+j]);
    states.push_back(state);
} // end function pushState

double Basis::jastrow(double a, double beta, double x12, double y12) {
    /* calculate Jastrow factor */
    return a/(beta + 1/sqrt(x12*x12 + y12*y12));
} // end function jastrow

double Basis::harmonicOscillatorWaveFunction(double alpha, double x, double y,
        int nx, int ny) {
    /* calculate harmonic oscillator wave function in 2D */
    return H(x,nx)*H(y,ny) * exp(-alpha*(x*x+y*y)/2);
} // end function harmonicOscillatorWaveFunction

void Basis::setBasisMatrix(Eigen::MatrixXd r, double alpha) {
    /* set matrix for Slater determinant with harmonic oscillator */
    double N = r.rows()/2;
    phiU.resize(N,N);
    phiD.resize(N,N);
    for (unsigned int i = 0; i < N; ++i) {
        for (unsigned int j = 0; j < N; ++j) {
            phiD(i,j) = harmonicOscillatorWaveFunction(alpha, r(2*i,0),
                    r(2*i,1), *states[2*j][0], *states[2*j][1]);
            phiU(i,j) = harmonicOscillatorWaveFunction(alpha, r(2*i+1,0),
                    r(2*i+1,1), *states[2*j+1][0], *states[2*j+1][1]);
        } // end forj
    } // end fori
} // end function setBasisMatrix

double Basis::trialWaveFunction(Eigen::MatrixXd r, double alpha, double beta,
        double a) {
    /* given a vector of coordinates, return trial wave function */
    unsigned int N = r.rows();
    double expInner = 0;
    setBasisMatrix(r,alpha);
    for (unsigned int i = 0; i < N; ++i) {
        for (unsigned int j = i+1; j < N; ++j) {
            expInner += jastrow(a, beta, r(i,0)-r(j,0),r(i,1)-r(j,1));
        } // end forj
    } // end fori
    return phiU.determinant() * phiD.determinant() * exp(expInner);
} // end function trialWaveFunction

void Basis::printStates() {
    /* print states */
    std::vector<int> nm(2);
    std::cout << "(nx,ny,s,ms,E,N)" << std::endl;
    for (unsigned int i = 0; i < states.size(); ++i) {
        std::cout << "(" 
            << *states[i][0] << ","
            << *states[i][1] << ","
            << *states[i][2] << ","
            << *states[i][3] << ","
            << *states[i][4] << ","
            << *states[i][5] << ")"
            << std::endl;
    } // end fori
    std::cout << "(nx,ny,s,ms,E,N)" << std::endl;
    std::cout << "Number of states: " << states.size() << std::endl;
} // end function print state

Basis::~Basis() {
    delete meth;
} // end deconstructor

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
#include <Eigen/Eigenvalues>

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

double Basis::jastrow(double x12, double y12) {
    /* calcualte Jastrow factor */
    return a/(beta + 1/sqrt(x12*x12 + y12*y12));
} // end function jastrow

double Basis::oscillator1D(double x, int n) {
    /* calculate hermite in 1D */
    return H(x,n) * exp(-alpha*x*x/2);
} // end function oscillator1D

double Basis::harmonicOscillatorWaveFunction(double x, double y, 
        int nx, int ny) {
    /* calculate harmonic oscillator wave function in 2D */
    return omega * oscillator1D(x,nx) * oscillator1D(y,ny);
} // end function harmonicOscillatorWaveFunction

double Basis::trialWaveFunction(Eigen::MatrixXd r, double alpha, double beta,
        double a) {
    /* given a vector of coordinates, return trial wave function */
    alpha = alpha;
    beta = beta;
    unsigned int N = r.rows();
    Eigen::MatrixXd Phi = Eigen::MatrixXd::Zero(N,N);
    double expInner = 0;
    double tmpWave = 0;
    for (unsigned int i = 0; i < N; ++i) {
        tmpWave = harmonicOscillatorWaveFunction(r(i,0), r(i,0), *states[i][0],
                *states[i][1]);
        for (unsigned int j = 0; j < N; ++j) {
            Phi(i,j) = tmpWave * harmonicOscillatorWaveFunction(r(j,0), r(j,0),
                    *states[j][0], *states[j][1]);
            if(i < j) {
                expInner += jastrow(r(i,0)-r(j,0),r(i,1)-r(j,1));
            } // end if
        } // end forj
    } // end fori
    return Phi.determinant() * exp(expInner);
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

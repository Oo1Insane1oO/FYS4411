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

double Basis::padejastrow(const unsigned int &i, const unsigned int &j) {
    /* return 1 for anti-parallel spin and 1/3 for parallel */
    if ((i%2 && j%2) || (!(i%2) && !(j%2))) {
        return 1./3;
    } else {
        return 1;
    } // end ifelse
} // end function padejastrow

double Basis::jastrow(const Eigen::MatrixXd &r, double beta) {
    /* calculate Jastrow factor */
    double factor = 0;
    for (unsigned int i = 0; i < r.rows(); ++i) {
        for (unsigned int j = i+1; j < r.rows(); ++j) {
            factor += padejastrow(i,j)/(beta + 1/sqrt(pow(r(i,0)-r(j,0),2) +
                        pow(r(i,1)-r(j,1),2)));
        } // end forj
    } // end fori
    return factor;
} // end function jastrow

double Basis::harmonicOscillatorWaveFunction(double alpha, double x, double y,
        int nx, int ny) {
    /* calculate harmonic oscillator wave function in 2D */
    return H(x,nx)*H(y,ny) * exp(-alpha*(x*x+y*y)/2);
} // end function harmonicOscillatorWaveFunction

void Basis::setTrialWaveFunction(Eigen::MatrixXd &psiD, Eigen::MatrixXd &psiU,
        const Eigen::MatrixXd &r, double alpha) {
    /* given a vector of coordinates, return trial wave function */
    for (unsigned int i = 0; i < r.rows(); i+=2) {
        for (unsigned int j = 0; j < r.rows(); j+=2) {
            psiD(i/2,j/2) = harmonicOscillatorWaveFunction(alpha, r(i,0), r(i,1),
                    *states[j][0], *states[j][1]);
            psiU(i/2,j/2) = harmonicOscillatorWaveFunction(alpha, r(i+1,0), r(i+1,1),
                    *states[j+1][0], *states[j+1][1]);
        } // end forj
    } // end fori
} // end function trialWaveFunction

void Basis::updateTrialWaveFunction(Eigen::MatrixXd &psi, const Eigen::MatrixXd
        &r, double alpha, unsigned int i) {
    /* update row i */
    for (unsigned int j = 0; j < r.rows()/2; ++j) {
        psi(i,j) = harmonicOscillatorWaveFunction(alpha, r(i,0), r(i,1),
                *states[j][0], *states[j][1]);
    } // end forj
} // end function updateTrialWaveFunction

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

std::vector<int> Basis::getMagicNumbers() {
    /* return vector containint magic numbers in acending order */
    return M;
} // end function getMagicNumbers

Basis::~Basis() {
    delete meth;
} // end deconstructor

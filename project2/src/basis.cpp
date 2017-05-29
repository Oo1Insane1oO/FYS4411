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
#include <numeric> // accumulate

Basis::Basis(double w, unsigned int cut) {
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

    // find magic numbers (given number of particles)
    int i = 1;
    M.insert(M.begin(),2);
    while (*(M.end()-1)<2*cut) {
        M.push_back(M[i-1] + 2*(i+1));
        i++;
    } // end while

    // set possible values for nx and ny and energy E
    n.resize(M.size());
    E.resize(M.size());
    for (unsigned int i = 0; i < n.size(); ++i) {
        n[i] = i;
        E[i] = i+1;
    } // end fori

    // allocate state arrays as {nx,ny,s,ms,E,M}
    std::vector<int*> s1 = std::vector<int*>(6,0);
    std::vector<std::vector<int>> tmpn;
    for (unsigned int i = 0; i < M.size(); ++i) {
        /* loop over states for energy level i */
        tmpn.clear();
        findPossiblenxny(i,tmpn);
        for (unsigned int j = 0; j < tmpn.size(); ++j) {
            pushState(s1, i, tmpn[j][0], tmpn[j][1], 0);
            pushState(s1, i, tmpn[j][0], tmpn[j][1], 1);
        } // end forj
    } // end fori
} // end constructor

void Basis::findPossiblenxny(unsigned int i, std::vector<std::vector<int>> &p) {
    /* find possible values of nx and ny given energy level */
    std::vector<int> nxny = std::vector<int>(2,0);
    for (unsigned int j = 0; j < E[i]; ++j) {
        for (unsigned int k = 0; k < E[i]; ++k) {
            if (j+k == E[i]-1) {
                nxny[0] = j;
                nxny[1] = k;
                p.push_back(nxny);
            } // end ifelse
        } // end fork
    } // end forj
} // end function findPossiblenxny

void Basis::pushState(std::vector<int*> &state, int e, int i, int j, int ud) {
    /* set (nx,ny,s,ms,E,M) in s1 and push s1 to states */
    state[0] = &(n[i]);
    state[1] = &(n[j]);
    state[2] = &s;
    state[3] = &(ms[ud]);
    state[4] = &(E[e]);
    state[5] = &(M[e]);
    states.push_back(state);
} // end function pushState

double Basis::padejastrow(const unsigned int &i, const unsigned int &j) {
    /* return 1 for anti-parallel spin and 1/3 for parallel */
    if (std::abs(i-j)>ECut) {
        return 1./3;
    } else {
        return 1;
    } // end if
} // end function padejastrow

double Basis::jastrow(const Eigen::MatrixXd &r, double beta) {
    /* calculate Jastrow factor */
    double factor = 0;
    for (unsigned int i = 0; i < r.rows(); ++i) {
        for (unsigned int j = i+1; j < r.rows(); ++j) {
            factor += padejastrow(i,j)/(beta + 1/(r.row(i)-r.row(j)).norm());
        } // end forj
    } // end fori
    return factor;
} // end function jastrow

void Basis::updateJastrow(double &factor, const Eigen::MatrixXd &rold, const
        Eigen::MatrixXd &rnew, double beta, unsigned int k) {
    /* update Jastrow factor for row k */
    for (unsigned int j = 0; j < rold.rows(); ++j) {
        if (j != k) {
            factor -= padejastrow(k,j)/(beta +
                    1/(rold.row(k)-rold.row(j)).norm());
            factor += padejastrow(k,j)/(beta +
                    1/(rnew.row(k)-rnew.row(j)).norm());
        } // end if
    } // end forj
} // end function updateJastrow

double Basis::jastrowRatio(const Eigen::MatrixXd &rold, const Eigen::MatrixXd
        &rnew, double beta, unsigned int k) {
    /* calculate ratio of jastrow factors (when only one row in Slater
     * determinant has changed) */
    double ratio = 0;
    for (unsigned int j = 0; j < rold.rows(); ++j) {
        if (j != k) {
            ratio += padejastrow(k,j) * (1/(beta +
                        1/(rnew.row(k)-rnew.row(j)).norm()) - 1/(beta +
                        1/(rold.row(k)-rold.row(j)).norm()));
        } // end if
    } // end forj
    return ratio;
} // end function jastrowRatio

double Basis::harmonicOscillatorWaveFunction(double alpha, double x, double y,
        int nx, int ny) {
    /* calculate harmonic oscillator wave function in 2D */
    return H(sqrt(omega*alpha)*x,nx) * H(sqrt(omega*alpha)*y,ny) *
        exp(-alpha*omega*(x*x+y*y)/2.0);
} // end function harmonicOscillatorWaveFunction

void Basis::setTrialWaveFunction(Eigen::MatrixXd &psiD, Eigen::MatrixXd &psiU,
        const Eigen::MatrixXd &r, const double alpha) {
    /* given a vector of coordinates, return trial wave function */
    for (unsigned int i = 0; i < r.rows()/2; ++i) {
        for (unsigned int j = 0; j < r.rows(); j+=2) {
            psiD(i,j/2) = harmonicOscillatorWaveFunction(alpha, r(i,0), r(i,1),
                    *states[j][0], *states[j][1]);
            psiU(i,j/2) = harmonicOscillatorWaveFunction(alpha,
                    r(i+r.rows()/2,0), r(i+r.rows()/2,1), *states[j+1][0],
                    *states[j+1][1]);
        } // end forj
    } // end fori
} // end function trialWaveFunction

void Basis::updateTrialWaveFunction(Eigen::MatrixXd &psi, const Eigen::MatrixXd
        &r, const double alpha, const unsigned int k, const unsigned int kIdx,
        const unsigned int jstart) {
    /* update row i given r(i) */
    for (unsigned int j = 0; j < r.rows(); j+=2) {
        psi(kIdx,j/2) = harmonicOscillatorWaveFunction(alpha, r(k,0), r(k,1),
                *states[j+jstart][0], *states[j+jstart][1]);
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
    /* return vector containing magic numbers in acending order */
    return M;
} // end function getMagicNumbers

Basis::~Basis() {
    delete meth;
} // end deconstructor

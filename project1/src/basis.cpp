//////////////////////////////////////////////////////////////////////////////
// Class for assembling a harmonic oscillator basis for calculating         //
// Hartree-Fock energies.                                                   //
//                                                                          //
// Main Functions:                                                          //
//      - HartreeFock: Hartree-Fock algortihm. Uses Eigen to find           // 
//                     eigenvalues and eigenvectors.                        //
//      - Assemble: Set interaction matrix elements. Uses an analytic       //
//                  expression to calculate the elements.                   //
//      - ConvertToPolar: Give principal quantum numbers in polar           // 
//                        coordinates given cartesian nx and ny.            //
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
#include "Coulomb_Functions.hpp"
#include <fstream>
#include <math.h>
#include <Eigen/Eigenvalues>

Basis::Basis(double w, int cut) {
    /* initialize states, set methods */
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

Basis::~Basis() {
    delete meth;
} // end deconstructor

unsigned int Basis::iidx(unsigned int N, unsigned int i, unsigned int j,
        unsigned int k, unsigned int l) {
    /* calculate index I[i,j,k,l] */
    return i + N * (j + N * (k + N*l));
} // end function index

void Basis::assemble(int t,bool perturb) {
    /* calculate interaction integrals and assemble antisymmetrized values */
    int pqS, pqM, rsM, rsS, tmpSizepqrs, tmpSizepqsr;
    std::array<int,2> MMSMap;
    if (t==0) {
        eps0Integrals.resize(states.size());
    } // end if
    std::vector<int> pstate, qstate, rstate, sstate;
    std::array<unsigned int,4> pqrs, pqsr, qpsr, qprs;
    double tmp;
    #pragma omp parallel for \
    private(tmp,MMSMap,pqS,pqM,rsM,rsS, \
            pstate,qstate,rstate,sstate, \
            pqrs, pqsr, qpsr, qprs, \
            tmpSizepqrs, tmpSizepqsr)
    for (unsigned int p = 0; p < states.size(); ++p) {
        if (t==0) {
            eps0Integrals[p] = omega * *states[p][4];
        } // end if
        pstate = convertToPolar(*states[p][0], *states[p][1]);
        pqrs[0] = p;
        pqsr[0] = p;
        qpsr[1] = p;
        qprs[1] = p;
        for (unsigned int q = p; q < states.size(); ++q) {
            qstate = convertToPolar(*states[q][0], *states[q][1]);
            pqM = pstate[1] + qstate[1];
            pqS = *states[p][3] + *states[q][3];
            MMSMap = {pqM,pqS};
            pqrs[1] = q;
            pqsr[1] = q;
            qpsr[0] = q;
            qprs[0] = q;
            for (unsigned int r = 0; r < states.size(); ++r) {
                rstate = convertToPolar(*states[r][0], *states[r][1]);
                pqrs[2] = r;
                pqsr[3] = r;
                qpsr[3] = r;
                qprs[2] = r;
                for (unsigned int s = r; s < states.size(); ++s) {
                    rsS = *states[r][3] + *states[s][3];
                    if (pqS != rsS) {
                        /* check spin conservation */
                        continue;
                    } // end if
                    sstate = convertToPolar(*states[s][0], *states[s][1]);
                    rsM = rstate[1] + sstate[1];
                    if (pqM != rsM) {
                        /* check orbital quantum number conservation */
                        continue;
                    } // end if
                    pqrs[3] = s;
                    pqsr[2] = s;
                    qpsr[2] = s;
                    qprs[3] = s;
                    # pragma omp critical 
                    {
                    if (t==0) {
                        /* set sizes and index-mappings */
                        tmpSizepqrs = sizes[MMSMap];
                        tmpSizepqsr = sizes[MMSMap] + 1;
                        Imap.emplace(std::make_pair(pqrs,tmpSizepqrs));
                        Imap.emplace(std::make_pair(pqsr,tmpSizepqsr));
                        Imap.emplace(std::make_pair(qpsr,tmpSizepqrs));
                        Imap.emplace(std::make_pair(qprs,tmpSizepqsr));
                        sizes[MMSMap] += 2;
                    } else if (t==1) {
                        /* resize vectors holding integral values */
                        integrals[MMSMap].resize(sizes[MMSMap],0);
                    } else {
                        /* calculate and set antisymetrized elements */
                        if (perturb) {
                            /* used for checking unperturb energies */
                            tmp = 0;
                        } else {
                            tmp = ((meth->kronk(*states[p][3], *states[r][3])
                                        && meth->kronk(*states[q][3],
                                            *states[s][3])) ?
                                    Coulomb_HO(omega, pstate[0], pstate[1],
                                        qstate[0], qstate[1], rstate[0],
                                        rstate[1], sstate[0], sstate[1]) : 0) -
                                ((meth->kronk(*states[p][3], *states[s][3]) &&
                                  meth->kronk(*states[q][3], *states[r][3])) ?
                                 Coulomb_HO(omega, pstate[0], pstate[1],
                                     qstate[0], qstate[1], sstate[0],
                                     sstate[1], rstate[0], rstate[1]) : 0);
                        } // end ifelse
                        integrals[MMSMap][Imap[pqrs]] = tmp;
                        integrals[MMSMap][Imap[pqsr]] = -tmp;
                    } // end if-elifelse
                    }
                } // end fors
            } // end forr
        } // end forq
    } // end forp
    if (t==0) {
        assemble(1,perturb);
    } else if (t==1) {
        assemble(2,perturb);
    } else {
        return;
    } // end ifelifelse
} // end function assemble

std::vector<int> Basis::convertToPolar(int nx, int ny) {
    /* map cartesian quantum numbers to polar */
    int nm[2] = {(nx>ny ? ny : nx), nx-ny};
    return std::vector<int> (std::begin(nm),std::end(nm));
} // end function convertToPolar

void Basis::setDensityMatrix(unsigned int n, Eigen::MatrixXd &densityMatrix,
        Eigen::MatrixXd C) {
    /* set density matrix in HartreeFock */
    for (unsigned int c = 0; c < C.rows(); ++c) {
        for (unsigned int d = 0; d < C.cols(); ++d) {
            densityMatrix(c,d) = 0;
            for (unsigned int i = 0; i < n; ++i) {
                densityMatrix(c,d) += C(c,i) * C(d,i);
            } // end fori
        } // end ford
    } // end forc
} // end function setDensityMatrix

void Basis::HartreeFock(unsigned int numParticles, int &maxIter, double eps) {
    /* Perform the Hartree-Fock algorithm and find ground state energy.
     * Interaction matrix needs to be assembled first.
     * Referenced maxIter is set to number of iterations needed for convergence
     * to be reached at the end*/
    double diff = 1. + eps;
    unsigned int N = states.size();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolve;
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(N,N);
    Eigen::MatrixXd densityMatrix = Eigen::MatrixXd::Zero(N,N);
    Eigen::VectorXd previousEnergies = Eigen::VectorXd::Zero(N);
    Eigen::MatrixXd HFmatrix = Eigen::MatrixXd::Zero(N,N);

    // initialize coefficients, assume diagonal and calculate density matrix
    setDensityMatrix(numParticles, densityMatrix,
            Eigen::MatrixXd::Identity(N,N));

    // run algorithm
    std::array<int,2> MMSMap;
    double sum = 0.0;
    int am, bm, cm, lS, lM;
    int start = 0;
    std::array<unsigned int,4> idx;
    while (start < maxIter && diff > eps) {
        /* perform HF-iteration */
        #pragma omp parallel for private(am,bm,lS,lM,MMSMap,idx) \
        reduction(+:sum)
        for (unsigned int a = 0; a < N; ++a) {
            am = convertToPolar(*states[a][0], *states[a][1])[1];
            idx[0] = a;
            for (unsigned int b = a; b < N; ++b) {
                bm = convertToPolar(*states[b][0], *states[b][1])[1];
                sum = 0;
                idx[2] = b;
                for (unsigned int c = 0; c < N; ++c) {
                    lS = *states[a][3] + *states[c][3];
                    lM = am + convertToPolar(*states[c][0],
                            *states[c][1])[1];
                    MMSMap = {lM,lS};
                    idx[1] = c;
                    for (unsigned int d = 0; d < N; ++d) {
                        if ((lS == *states[b][3]+*states[d][3]) && (lM == bm +
                                    convertToPolar(*states[d][0],
                                        *states[d][1])[1])) {
                            /* check orbital- and spin quantum number
                             * conservation */
                            idx[3] = d;
                            sum += densityMatrix(c,d) *
                                integrals[MMSMap][Imap[idx]];
                        } // end if
                    } // end ford
                } // end forc

                // add single particle energy term to diagonal
                HFmatrix(a,b) = ((a==b) ? (sum + eps0Integrals[a]) : sum);

                // HFmatrix is symmetric by definition
                HFmatrix(b,a) = HFmatrix(a,b);
            } // end forb
        } // end fora

        // compute eigenvalues(energies), eigenvectors(coefficients) and new
        // density matrix
        eigenSolve.compute(HFmatrix);
        C = eigenSolve.eigenvectors();
        setDensityMatrix(numParticles, densityMatrix,
                eigenSolve.eigenvectors());

        // set value for convergence test
        diff = 0;
        for (unsigned int i = 0; i < N; ++i) {
            diff += std::fabs(eigenSolve.eigenvalues()(i) -
                    previousEnergies(i));
        } // end fori
        diff /= N;

        // update energies, increment HF count
        previousEnergies = eigenSolve.eigenvalues();
        start++;
    } // end while

    // set single particle energies in class
    singleParticleEnergiesHartreeFock.resize(eigenSolve.eigenvalues().size());
    for (unsigned int i = 0; i < eigenSolve.eigenvalues().size(); ++i) {
        singleParticleEnergiesHartreeFock[i] = eigenSolve.eigenvalues()(i);
    } // end fori

    // find ground state energy estimation
    double E0 = 0;
    for (unsigned int i = 0; i < numParticles; ++i) {
        E0 += eigenSolve.eigenvalues()(i);
    } // end fori

    # pragma omp parallel for private(am,idx,lS,lM,MMSMap,cm) reduction(-:E0)
    for (unsigned int a = 0; a < N; ++a) {
        am = convertToPolar(*states[a][0], *states[a][1])[1];
        idx[0] = a;
        for (unsigned int b = 0; b < N; ++b) {
            lS = *states[a][3] + *states[b][3];
            lM = am + convertToPolar(*states[b][0],
                    *states[b][1])[1];
            MMSMap = {lM,lS};
            idx[1] = b;
            for (unsigned int c = 0; c < N; ++c) {
                cm = convertToPolar(*states[c][0], *states[c][1])[1];
                idx[2] = c;
                for (unsigned int d = 0; d < N; ++d) {
                    if ((lS == *states[c][3]+*states[d][3]) && (lM ==
                                cm+convertToPolar(*states[d][0],
                                    *states[d][1])[1])) { 
                        idx[3] = d;
                        E0 -= 0.5 * densityMatrix(a,c) * densityMatrix(b,d) *
                            integrals[MMSMap][Imap[idx]];
                    } // end if
                } // end ford
            } // end forc
        } // end forb
    } // end fora
    E0HartreeFock = E0;
    maxIter = start;
} // end function HartreeFock

void Basis::printStates() {
    /* print states */
    std::vector<int> nm(2);
    std::cout << "(nx,ny,s,ms,E,N) (n,m)" << std::endl;
    for (unsigned int i = 0; i < states.size(); ++i) {
        nm = convertToPolar(*states[i][0], *states[i][1]);
        std::cout << "(" 
            << *states[i][0] << ","
            << *states[i][1] << ","
            << *states[i][2] << ","
            << *states[i][3] << ","
            << *states[i][4] << ","
            << *states[i][5] << ")"
            << " (" << nm[0] << "," << nm[1] << ")" << " " << i
            << std::endl;
    } // end fori
    std::cout << "(nx,ny,s,ms,E,N)" << std::endl;
    std::cout << "Number of states: " << states.size() << std::endl;
} // end function print state

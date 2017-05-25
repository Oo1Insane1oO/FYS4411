//////////////////////////////////////////////////////////////////////////////
// Class containing functions for testinc class Basis.						//
//																			//
// All functions retunr true in case where the tests are same and false		//
// otherwise. Except for run_test which just runs all the functions.		//
//////////////////////////////////////////////////////////////////////////////

#include "tests.h" // header
#include <iostream>
#include <iomanip>

Tests::Tests(Basis *B, VMC *V, int n) {
    b = B;
    v = V;
    m = new Methods();

    eps = 1e-14;

    std::mt19937_64 mt(123);
    std::uniform_real_distribution<double> dist(0,1);

    oldM = Eigen::MatrixXd::Zero(n,n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            oldM(i,j) = dist(mt);
        } // end forj
    } // end fori
    rowi = floor(static_cast<double>(n)/2.);
    newM = oldM;
    for (int j = 0; j < n; ++j) {
        newM(rowi,j) = dist(mt);
    } // end forj
} // end constructor

Tests::~Tests() {
    delete b;
    delete v;
    delete m;
} // end deconstructor

bool Tests::test_energies() {
    /* check that energy in state is correct */
    bool t = false;
    for (unsigned int i = 0; i < b->states.size(); ++i) {
        if (*(b->states)[i][4] == *(b->states)[i][0]+*(b->states)[i][1] + 1) {
            t = true;
        } else {
            t = false;
        } // emd ifelse
        if (!t) {
            break;
        } // end if
    } // end fori
    return t;
} // end function test_energies

bool Tests::test_2particle() {
    /* check that energy in case of unperturbed harmonic oscillator system with
     * 2 electrons is correct */
    v->setCoulombInteraction(false);
    v->calculate();
    return ((fabs(v->energy-2)<=eps &&
                fabs(m->variance(v->energy,v->energySq))<=eps) ? true : false);
} // end function test_2particle

bool Tests::test_determinantratio() {
    /* check that ratio between determinants are correct in member function
     * determinantRatio of class methods */
    return (std::fabs(m->determinantRatio(newM,oldM.inverse(),rowi) -
                newM.determinant()/oldM.determinant())<=1e-13 ? true : false);
} // end function test_determinantratio

bool Tests::test_updateinverse() {
    /* check that inverse is updated correctly */
    Eigen::MatrixXd newInv = Eigen::MatrixXd::Zero(oldM.rows(),oldM.rows());
    Eigen::MatrixXd inv = newM.inverse();
    m->updateMatrixInverse(oldM, newM, oldM.inverse(), newInv,
            m->determinantRatio(newM, oldM.inverse(), rowi), rowi);
    bool t = false;
    for (int i = 0; i < oldM.rows(); ++i) {
        for (int j = 0; j < oldM.rows(); ++j) {
            if (std::fabs(newInv(i,j) - inv(i,j))<=1e-13) {
                t = true;
            } else {
                t = false;
                break;
            } // end if
        } // end forj
        if (!t) {
            break;
        } // end if
    } // end fori
    return t;
} // end function test_updateinverse

bool Tests::test_wavefunction2() {
    /* test wavefunction for 2 electrons */
    Eigen::MatrixXd oldD = Eigen::MatrixXd::Zero(b->ECut,b->ECut);
    Eigen::MatrixXd oldU = Eigen::MatrixXd::Zero(b->ECut,b->ECut);
    b->setTrialWaveFunction(oldD, oldU, oldM, v->alpha);
    double trial = oldD.determinant() * oldU.determinant();
    double wave2 = b->harmonicOscillatorWaveFunction(v->alpha, oldM(0,0),
            oldM(0,1), 0, 0) * b->harmonicOscillatorWaveFunction(v->alpha,
            oldM(1,0), oldM(1,1), 0, 0);
    return (fabs(trial-wave2) <= eps ? true : false);
} // end function test_wavefunction2

bool Tests::test_padejastrow() {
    /* check that Pade-Jastrow factor gives correct values for corresponding
     * odd and even states */
    bool t = false;
    double a;
    for (unsigned int i = 0; i < b->states.size(); ++i) {
        for (unsigned int j = i+1; j < b->states.size(); ++j) {
            a = b->padejastrow(i,j);
            if (((i%2 && j%2) || (!(i%2) && !(j%2))) && (a-1/3.)<=eps) {
                t = true;
            } else if (((i%2 && !(j%2)) || (!(i%2) && j%2)) && (a-1)<=eps) {
                t = true;
            } else {
                t = false;
                break;
            } // end ifelse
        } // end forj
        if (!t) {
            break;
        } // end if
    } // end foi
    return t;
} // end function test_padejastrow

bool Tests::test_conjugateGradient() {
    /* check ouput from conjugateGradient in class Methods */
    // set A and b in Ax=b matrix
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2,2);
    Eigen::MatrixXd startx = Eigen::MatrixXd::Zero(2,1);
    Eigen::MatrixXd b = Eigen::MatrixXd::Zero(2,1);
    A << 4, 1,
         1, 3;
    startx << 2, 1;
    b << 1, 2;
    Eigen::MatrixXd x = m->conjugateGradient(A, b, startx);
    return ((fabs(x(0)-1./11)<=eps && fabs(x(1)-7./11)<=eps ? true : false));
} // end function test_conjugateGradient

void Tests::run_tests(int t) {
    /* run all tests and exit */
    if (t) {
        if(test_energies()) {
            std::cout << "Energies good" << std::endl;
        } else {
            std::cout << "Energies wrong" << std::endl;
        } // end ifelse
        if(test_2particle()) {
            std::cout << "Energy unperturbed 2 electron good" << std::endl;
            std::cout << std::setprecision(10) << "  Energy is: " << v->energy << std::endl;
        } else { 
            std::cout << "Energy unperturbed 2 electron wrong" << std::endl;
            std::cout << std::setprecision(10) << "  Energy is: " << v->energy << std::endl;
        } // end ifelse
        if(test_determinantratio()) {
            std::cout << "Determinant ratio good" << std::endl;
        } else { 
            std::cout << "Determinant ratio wrong" << std::endl;
        } // end ifelse
        if(test_updateinverse()) {
            std::cout << "Update inverse good" << std::endl;
        } else { 
            std::cout << "Update inverse wrong" << std::endl;
        } // end ifelse
        if (b->ECut == 1) {
            if(test_wavefunction2()) {
                std::cout << "Wavefunction 2 electron good" << std::endl;
            } else { 
                std::cout << "Wavefunction 2 electron wrong" << std::endl;
            } // end ifelse
        } // end if
        if(test_padejastrow()) {
            std::cout << "Pade-Jastrow factor good" << std::endl;
        } else {
            std::cout << "Pade-Jastrow factor wrong" << std::endl;
        } // end ifelse
        if(test_conjugateGradient()) {
            std::cout << "Conjugate Gradient good" << std::endl;
        } else {
            std::cout << "Conjugate Gradient wrong" << std::endl;
        } // end ifelse
        if (t==2) {
            b->printStates();
        } // end ift
        exit(1);
    } // end if t
} // end function run_tests

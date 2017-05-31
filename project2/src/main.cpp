#include "tests.h" // test functions#include "vmc.h" // class basis
#include <stdlib.h> // atoi
#include <iostream> // cout
#include <chrono> // timer
#include <iomanip> // setprecision
#include <algorithm> // find
#include <chrono> // timer
#include <sstream>
#include <iterator>
#include <string.h>
#include <mpi.h>

//////////////////////////////////////////////////////////////////////////////
// Main file for running vmc algorithm                                      //
//////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
    /* main */

    int myRank, numProcs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Request req;

    if (argc < 9) {
        /* Print usage if number of command line arguments are to few */
        std::cout << 
            "USAGE: mpirun -np <number of processes> ./main 'omega' 'particles' 'iterations' 'tests' 'importance' 'Coulomb' 'Jastrow' 'num procs" 
            << std::endl;
        std::cout <<
            "    " << "omega: (float) HO frequency\n" <<
            "    " << "particles: (int) Fermi level(closed shell)\n" <<
            "    " << "iterations: (int) Max iterations in VMC(MC cycles) \n" <<
            "    " << "step: (float) step size in VMC \n" <<
            "    " << "tests: (1/0) indicating to run tests or not\n" <<
            "    " << "importance sampling: (1/0) indicating to run with importance sampling or not\n" <<
            "    " << "Coulomb: (1/0) Indicating to run with Coulomb interaction or not\n" <<
            "    " << "Jastrow: (1/0) Indicating to run with Jastrow factor or not" <<
            "    " << "Destinaiton: (string) Destination folder for saving files (can be ignored)" <<
            std::endl;
        exit(1);
    } else if (numProcs < 1) {
        std::cout << "Number of processes needs to be at least 1" << std::endl;
        MPI_Finalize();
        exit(1);
    } // end if

    // grab parameters as command line arguments
    double omega = atof(argv[1]);
    int num = atoi(argv[2]);
    unsigned int maxIterations = atoi(argv[3]);
    double step = atof(argv[4]);
    int t = atoi(argv[5]);
    bool imp = atoi(argv[6]);
    bool coul = atoi(argv[7]);
    bool jast = atoi(argv[8]);
    const char *filename;
    if (argc == 10) {
        filename = argv[9];
    } else {
        filename = NULL;
    } // end ifelse

    // let Eigen use openmp
    Eigen::initParallel();

    // divide number of variational runs between the processes evenly
    unsigned int maxCount = 1000;
    float tmpNum = (float)maxCount / numProcs;
    unsigned int myMaxCount = (myRank < maxCount % numProcs ? ceil(tmpNum) :
            floor(tmpNum));
    
    // set basis (cartesian)
    Basis *b = new Basis(omega, num/2);

    // make sure number of particles is a magic number(closed shell)
    std::vector<int> magicNumber = b->getMagicNumbers();
    std::vector<int>::iterator it;
    it = std::find(magicNumber.begin(), magicNumber.end(), num);
    if (it == magicNumber.end()) {
        std::cout << "make sure num is a magic number N=2,6,12,20,30,42..." <<
            std::endl;
        exit(1);
    } // end if
    
    // set vmc object for calculations (with different seed for each process)
    double myAlpha;
    double myBeta;
    long long int mySeed;
    long long int seedbuf;
    double alphaBuf;
    double betaBuf; 
    if (myRank == 0) {
        /* let master distribute seeds and initial variational parameters(and
         * make his own) */
        std::istringstream stringBuffer("0 1 2 3 4 5 6 7 8 9 10");
        std::istream_iterator<int> start(stringBuffer), end;
        std::seed_seq seedSequence(start, end);
        std::mt19937_64 generator(seedSequence);
        std::uniform_real_distribution<double> dist(0.1,1.3);
        mySeed = std::chrono::high_resolution_clock::now() .
            time_since_epoch().count();
        myAlpha = dist(generator);
        myBeta = dist(generator);
        for (int p = 1; p < numProcs; ++p) {
            seedbuf = std::chrono::high_resolution_clock::now() .
                time_since_epoch().count();
            alphaBuf = dist(generator);
            betaBuf = dist(generator);
            MPI_Send(&seedbuf, 1, MPI_LONG_LONG, p, 0, MPI_COMM_WORLD);
            MPI_Send(&alphaBuf, 1, MPI_DOUBLE, p, 1, MPI_COMM_WORLD);
            MPI_Send(&betaBuf, 1, MPI_DOUBLE, p, 2, MPI_COMM_WORLD);
        } // end forp
    } else {
        /* Slaves receive seed */
        MPI_Recv(&mySeed, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
        MPI_Recv(&myAlpha, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
        MPI_Recv(&myBeta, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    } // end ifelse

    VMC *vmcObj = new VMC(b, myAlpha, myBeta, 2, step, maxIterations, mySeed);
    vmcObj->setImportanceSampling(imp);
    vmcObj->setCoulombInteraction(coul);
    vmcObj->setJastrow(jast);
    
    if (t) {
        /* run tests */
        Tests testObj = Tests(b,vmcObj,num);
        testObj.run_tests(t);
        exit(1);
    } // end if

    // run calculations
//     std::chrono::steady_clock::time_point begin;
//     begin = std::chrono::steady_clock::now();

    std::cout << "Starting :" << myRank << " " << myAlpha << " " << myBeta <<
        std::endl;
    // create filename for each process and run
    char myFileName[80];
    sprintf(myFileName, "%sP%d_", filename, myRank);
    vmcObj->calculate(myMaxCount, myFileName);

//     std::chrono::steady_clock::time_point end;
//     end = std::chrono::steady_clock::now();
//     std::cout << "Calculation time: " <<
//         std::chrono::duration_cast<std::chrono::seconds>(end-begin).count()
//         << std::endl;

//     std::cout << std::setprecision(10) << "<E> = " << vmcObj->energy << ", " <<
//         "<E^2> = " << vmcObj->energySq << std::endl;
//     std::cout << std::setprecision(10) << "<E^2> - <E>^2 = " <<
//         (vmcObj->energySq - pow(vmcObj->energy,2))/maxIterations << std::endl;
// 
//     std::cout << "alpha: " << vmcObj->alpha << ", beta: " << vmcObj->beta <<
//         std::endl;

    // free objects
    delete b;
    delete vmcObj;
    if (!filename) {
        delete filename;
    } // end if

    MPI_Finalize();

    return 0;
} // end main

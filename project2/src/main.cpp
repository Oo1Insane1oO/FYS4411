#include "tests.h" // test functions
#include "basis.h" // class Basis
#include "vmc.h" // class VMC
#include <stdlib.h> // atoi
#include <iostream> // cout
#include <chrono> // timer
#include <iomanip> // setprecision
#include <algorithm> // find
#include <chrono> // timer
#include <sstream> // for seeding
#include <iterator> // for seeding
#include <string.h> // for filename
#include <mpi.h> // header for MPI

//////////////////////////////////////////////////////////////////////////////
// Main file for running vmc algorithm                                      //
//                                                                          //
// Also sets up MPI and divides number of VMC cycles evenly between the     //
// processes and distributes a random seed to each process.                 //
//////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
    /* main */
    int myRank, numProcs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

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
            "    " << "Jastrow: (1/0) Indicating to run with Jastrow factor or not\n" <<
            "    " << "Destination: (string) Destination folder and filename for saving files (will not write to file if ignored)" <<
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
    if (argc >= 10) {
        filename = argv[9];
    } else {
        filename = NULL;
    } // end ifelse

    // set filename to NULL if empty
    if (filename) {
        std::string tmpf;
        tmpf = filename;
        if (!tmpf.compare(" ")) {
            filename = NULL;
        } // end if
    } // end if

    // let Eigen use openmp
    Eigen::initParallel();

    // divide number of variational runs between the processes evenly
    int maxCount = 200;
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
        MPI_Finalize();
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
//         std::istringstream stringBuffer("0 1 2 3 4 5 6 7 8 9 10");
//         std::istream_iterator<int> start(stringBuffer), end;
//         std::seed_seq seedSequence(start, end);
        mySeed = std::chrono::high_resolution_clock::now() .
            time_since_epoch().count() * myRank;
        std::mt19937_64 generator(mySeed);
        std::uniform_real_distribution<double> dista(0.8,1.0);
        std::uniform_real_distribution<double> distb(0.25,0.45);
//         std::uniform_real_distribution<double> dista(1.0,1.07);
//         std::uniform_real_distribution<double> distb(0.469,0.475);
        mySeed = std::chrono::high_resolution_clock::now() .
            time_since_epoch().count();
        myAlpha = dista(generator);
        myBeta = distb(generator);
        for (int p = 1; p < numProcs; ++p) {
            seedbuf = std::chrono::high_resolution_clock::now() .
                time_since_epoch().count();
            alphaBuf = dista(generator);
            betaBuf = distb(generator);
            MPI_Send(&seedbuf, 1, MPI_LONG_LONG, p, 0, MPI_COMM_WORLD);
            MPI_Send(&alphaBuf, 1, MPI_DOUBLE, p, 1, MPI_COMM_WORLD);
            MPI_Send(&betaBuf, 1, MPI_DOUBLE, p, 2, MPI_COMM_WORLD);
        } // end forp
    } else {
        /* Slaves receive */
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

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(1);
    } // end if

    // Run Monte Carlo simulations and find optimal parameters
    std::chrono::steady_clock::time_point begin;
    begin = std::chrono::steady_clock::now();
    if (argc < 11) {
        vmcObj->calculate(myMaxCount);
    } // end if

    // root sends new parameters(averaged) to all
    double newAlpha = 0;
    double newBeta = 0;
    MPI_Reduce(&(vmcObj->alpha), &newAlpha, 1, MPI_DOUBLE, MPI_SUM, 0,
            MPI_COMM_WORLD);
    MPI_Reduce(&(vmcObj->beta), &newBeta, 1, MPI_DOUBLE, MPI_SUM, 0,
            MPI_COMM_WORLD);
    if (myRank == 0) {
        newAlpha /= numProcs;
        newBeta /= numProcs;
    } // end if
    MPI_Bcast(&newAlpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&newBeta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // dont minimize, take parameters as input
    if(argc == 12) {
        newAlpha = atof(argv[10]);
        newBeta = atof(argv[11]);
    } // end if

    // set optimal parameters
    if (myRank == 0) {
        std::cout << std::setprecision(16) << "Alpha: " << newAlpha << " Beta: \
            "<< newBeta << std::endl;
    } // end if
    vmcObj->setAlpha(newAlpha);
    vmcObj->setBeta(newBeta);

    // divide number of samples evenly among processes
    maxIterations *= 100;
    double tmpMaxNum = (double)maxIterations / numProcs;
    unsigned int myMaxIteration = ((myRank < (maxIterations%numProcs)) ?
            std::ceil(tmpMaxNum) : std::floor(tmpMaxNum));
    vmcObj->maxIterations = myMaxIteration;
    
    // create filename for each process
    char myFileName[100];
    if (filename) {
        sprintf(myFileName, "%sP%d", filename, myRank);
    } // end fi

    // run last simulation and write to file
    vmcObj->calculate(1, myFileName);

    std::chrono::steady_clock::time_point end;
    end = std::chrono::steady_clock::now();
    double myTime = std::chrono::duration_cast<std::chrono::milliseconds>(end -
            begin).count();

    // gather and average energies from all processes
    double energy, energySq, acceptance;
    MPI_Reduce(&(vmcObj->energy), &energy, 1, MPI_DOUBLE, MPI_SUM, 0,
            MPI_COMM_WORLD);
    MPI_Reduce(&(vmcObj->energySq), &energySq, 1, MPI_DOUBLE, MPI_SUM, 0,
            MPI_COMM_WORLD);
    MPI_Reduce(&(vmcObj->acceptance), &acceptance, 1, MPI_DOUBLE, MPI_SUM, 0,
            MPI_COMM_WORLD);
    double totalTime;
    MPI_Reduce(&(myTime), &totalTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myRank == 0) {
        energy /= numProcs;
        energySq /= numProcs;
        totalTime /= numProcs;
        acceptance /= numProcs;

        if (totalTime >= 100) {
            std::cout << "Calculation time: " << totalTime*0.001 << "s" <<
                std::endl;
        } else if (totalTime >= 6e4) {
            std::cout << "Calculation time: " << totalTime*0.001/60 << "min" <<
                std::endl;
        } else if (totalTime >= 3.6e6) {
            std::cout << "Calculation time: " << totalTime*0.001/3600 << "h" <<
                std::endl;
        } else {
            std::cout << "Calculation time: " << totalTime << "ms" <<
                std::endl;
        } // end ifeifeifelse

        std::cout << std::setprecision(16) << "Acceptance: " << acceptance <<
            std::endl;
        std::cout << std::setprecision(16) << "<E> = " << energy << ", " <<
            "<E^2> = " << energySq << std::endl;
        std::cout << std::setprecision(16) << "var(E) = " << (energySq -
                pow(energy,2))/maxIterations << std::endl;

        std::cout << std::setprecision(16) << "alpha: " << vmcObj->alpha << ",\
            beta: " << vmcObj->beta << std::endl;
    } // end if

    // free objects
    delete b;
    delete vmcObj;

    // clean up MPI
    MPI_Finalize();

    return 0;
} // end main

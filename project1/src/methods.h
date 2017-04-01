//////////////////////////////////////////////////////////////////////////////
// Header file for class Methods, see .cpp file for more information        //
//////////////////////////////////////////////////////////////////////////////

#ifndef METHODS_H
#define METHODS_H

#include <vector>
#include <functional>

class Methods {
    private:
        void copyVector(std::vector<double>&, std::vector<double>);
        void setPoints(int, std::vector<double>&);
        void setWeight(int, std::vector<double>&, std::vector<double>&);
        double hermiteNormal(int);
        int checkn(int);

    public:
        Methods();

        int kronk(int*, int*);
        int kronk(int&, int&);
        int factorial(int);
        double hermite(double, int);
        double gaussHermiteQuadrature(int*,
                std::function<double(double,double,double,double)>);
        void findRootsHermite(int, std::vector<double>&);
        void generateRandom(std::vector<double>&);
};

#endif /* METHODS_H */

//////////////////////////////////////////////////////////////////////////////
// Header file for class Basis, see .cpp file for more information            //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#ifndef BASIS_H
#define BASIS_H

#include <Eigen/Dense>
#include <unordered_map>
#include <boost/unordered_map.hpp>
#include <array>

class Basis {
    private:
        int ECut;
        double omega;

    public:
        Basis(double, int);
        virtual ~Basis();
        
        unsigned int iidx(unsigned int, unsigned int, unsigned int, unsigned
                int, unsigned int);

        int s;
        std::vector<int> n, ms, E, M, m;
        std::vector<std::vector<int*>> states;
};

#endif /* BASIS_H */

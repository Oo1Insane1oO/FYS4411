//////////////////////////////////////////////////////////////////////////////
// Header file for class Basis, see .cpp file for more information          //
//                                                                          //
// Includes a template for creating an unordered_map (hash-table)           //
//////////////////////////////////////////////////////////////////////////////

#ifndef BASIS_H
#define BASIS_H

#include "methods.h"
#include <Eigen/Dense>
#include <unordered_map>
#include <boost/unordered_map.hpp>
#include <array>

class Basis {
    private:
        // template hash used in std::unordered_map
        template <typename seq> struct seq_hash {
            std::size_t operator() (const seq& s) const {
                // use boost to grab range of hash
                std::size_t hash = 0;
                boost::hash_range(hash,s.begin(),s.end());
                return hash;
            }
        };

        // template defining std::unordered_map using above template as hash
        template <typename seq, typename T>
        using stdumap = std::unordered_map<seq,T,seq_hash<seq>>;

        int ECut;
        double omega;

        Methods *meth;
        void setDensityMatrix(unsigned int, Eigen::MatrixXd&, Eigen::MatrixXd);

    public:
        Basis(double, int);
        virtual ~Basis();
        
        unsigned int iidx(unsigned int, unsigned int, unsigned int, unsigned
                int, unsigned int);

        int s;
        double E0HartreeFock;
        std::vector<double> singleParticleEnergiesHartreeFock;
        std::vector<int> n, ms, E, M, m;
        std::vector<std::vector<int*>> states;
        stdumap<std::array<int,2>, std::vector<double>> integrals;
        stdumap<std::array<unsigned int,4>, unsigned int> Imap;
        stdumap<std::array<int,2>, unsigned int> sizes;
        std::vector<double> eps0Integrals;

        void assemble(int=0,bool=false);
        void printStates();
        void HartreeFock(unsigned int, int&, double);
        std::vector<int> convertToPolar(int, int);
};

#endif /* BASIS_H */

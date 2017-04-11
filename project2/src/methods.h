//////////////////////////////////////////////////////////////////////////////
// Header file for class Methods, see .cpp file for more information        //
//////////////////////////////////////////////////////////////////////////////

#ifndef METHODS_H
#define METHODS_H

#include <vector>
#include <functional>
#include <Eigen/Dense>

class Methods {
    private:

    public:
        Methods();

        int factorial(int);
        void updateMatrixInverse(Eigen::MatrixXd, Eigen::MatrixXd,
                Eigen::MatrixXd, Eigen::MatrixXd&, unsigned int);
        double determinantRatio(Eigen::MatrixXd, Eigen::MatrixXd, unsigned
                int);
};

#endif /* METHODS_H */

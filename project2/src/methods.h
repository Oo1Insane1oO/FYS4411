//////////////////////////////////////////////////////////////////////////////
// Header file for class Methods, see .cpp file for more information        //
//////////////////////////////////////////////////////////////////////////////

#ifndef METHODS_H
#define METHODS_H

#include <vector>
#include <functional>

class Methods {
    private:

    public:
        Methods();

        int factorial(int);
        double updateMatrixInverse(Eigen::MatrixXd, Eigen::MatrixXd,
                Eigen::MatrixXd&);
        double determinantRatio(Eigen::MatrixXd, Eigen::MatrixXd);
};

#endif /* METHODS_H */

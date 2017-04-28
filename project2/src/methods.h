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
        double min(double, double);
        int min(int, int);
        double max(double, double);
        int max(int, int);
        void updateMatrixInverse(const Eigen::MatrixXd&, const
                Eigen::MatrixXd&, const Eigen::MatrixXd&, Eigen::MatrixXd&,
                unsigned int);
        double determinantRatio(const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                unsigned int);
        double variance(double, double);

        Eigen::MatrixXd conjugateGradient(Eigen::MatrixXd, Eigen::MatrixXd,
                Eigen::MatrixXd);
};

#endif /* METHODS_H */

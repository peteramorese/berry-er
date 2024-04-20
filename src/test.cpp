#include "Dynamics.h"

#include <iostream>

#include <Eigen/Dense>

using namespace BRY;

int main() {

    PolynomialDynamics<2> dynamics(1, 1);    
    dynamics[0].coeff(0, 0) = 5;
    dynamics[0].coeff(1, 0) = 3;
    dynamics[0].coeff(1, 1) = 2;

    dynamics[1].coeff(0, 0) = 6;
    dynamics[1].coeff(1, 0) = 4;
    dynamics[1].coeff(1, 1) = 7;

    DEBUG("Dynamics[0]: " << dynamics[0]);
    DEBUG("Dynamics[1]: " << dynamics[1]);

    Eigen::MatrixXd dyn_mat = dynamics.dynamicsPowerMatrix(2);
    Eigen::MatrixXd nice_dyn_mat = dyn_mat.unaryExpr([](double x) {
        if (std::abs(x) < 0.0001)
            return 0.0;
        else
            return x;
    });
    DEBUG("dyn mat: \n" << nice_dyn_mat);
    DEBUG("col: " << dynamics[0]);
    DEBUG("col: " << dynamics[0]);
    DEBUG("col: " << (dynamics[0]^2));
    DEBUG("col: " << (dynamics[0]^2) * (dynamics[1]^2));
    //for (auto midx = mIdx(2, 3); !midx.last(); ++midx) {
    //    INFO("midx: " << midx);
    //}
    return 0;
}
 
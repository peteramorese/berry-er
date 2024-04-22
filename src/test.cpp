#include "Dynamics.h"
#include "Noise.h"
#include "Transform.h"
#include "HyperRectangle.h"

#include <iostream>

#include <Eigen/Dense>

using namespace BRY;

int main() {

    auto printVec = [](const auto& vec) {
        for (std::size_t i = 0; i < vec.size(); ++i) {
            std::cout << vec(i);
            if (i < vec.size() - 1)
                std::cout << ", ";
        }
        std::cout << std::endl;
    };

    //PolynomialDynamics<2> dynamics(1, 1);    
    //dynamics[0].coeff(0, 0) = 5;
    //dynamics[0].coeff(1, 0) = 3;
    //dynamics[0].coeff(1, 1) = 2;

    //dynamics[1].coeff(0, 0) = 6;
    //dynamics[1].coeff(1, 0) = 4;
    //dynamics[1].coeff(1, 1) = 7;

    //DEBUG("Dynamics[0]: " << dynamics[0]);
    //DEBUG("Dynamics[1]: " << dynamics[1]);

    //Eigen::MatrixXd dyn_mat = dynamics.dynamicsPowerMatrix(2);
    //Eigen::MatrixXd nice_dyn_mat = dyn_mat.unaryExpr([](double x) {
    //    if (std::abs(x) < 0.0001)
    //        return 0.0;
    //    else
    //        return x;
    //});
    //DEBUG("dyn mat: \n" << nice_dyn_mat);
    //DEBUG("col: " << dynamics[0]);
    //DEBUG("col: " << dynamics[0]);
    //DEBUG("col: " << (dynamics[0]^2));
    //DEBUG("col: " << (dynamics[0]^2) * (dynamics[1]^2));

    //for (auto midx = mIdx(2, 3); !midx.last(); ++midx) {
    //    INFO("midx: " << midx);
    //}

    //Covariance<2> cov;
    //cov << 2.2, 1.1, 
    //       1.1, 3.3;
    
    //Additive2ndMomentNoise<2> noise(cov);

    //Eigen::MatrixXd add_noise_mat = noise.getAdditiveNoiseMatrix(2);
    //DEBUG("noise mat: \n" << add_noise_mat);

    //HyperRectangle<1> set;
    //set.lower_bounds = Eigen::Vector<bry_float_t, 1>{0.4};
    //set.upper_bounds = Eigen::Vector<bry_float_t, 1>{0.6};

    //Polynomial<1> p(4);

    //p.coeff(0) = 1;
    //p.coeff(1) = -1.7;
    //p.coeff(2) = 0.8;
    //p.coeff(3) = 0.0;
    //p.coeff(4) = 2.0;

    HyperRectangle<2> set;
    set.lower_bounds = {0.5, 0.5};
    set.upper_bounds = {1.0, 1.0};

    Polynomial<2> p(1);

    p.coeff(0, 0) = 1;
    p.coeff(1, 0) = 2;
    p.coeff(0, 1) = 3;
    p.coeff(1, 1) = -2;

    Eigen::MatrixXd tf = Transform::affine(set, p.degree());

    DEBUG("transform mat: \n" << tf);

    DEBUG("p before: " << p);
    //std::cout << p.tensor().transpose() << std::endl;
    printVec(p.tensor());
    auto p_tf = transform(p, tf);
    DEBUG("p after: " << p_tf);
    printVec(p_tf.tensor());
    //std::cout << p_tf.tensor().transpose() << std::endl;

    return 0;
}
 
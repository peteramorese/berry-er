#include "Dynamics.h"
#include "Noise.h"
#include "HyperRectangle.h"
#include "Synthesis.h"
#include "ArgParser.h"

#include <iostream>
#include <stdio.h>

#include <Eigen/Dense>

using namespace BRY;

int main(int argc, char** argv) {



    //Covariance<2> cov;
    //cov(0, 0) = 1.00;
    //cov(1, 0) = 6.00;
    //cov(0, 1) = 6.00;
    //cov(1, 1) = 4.00;

    //Eigen::Vector<bry_float_t, 2> mean;
    //mean(0) = 10;
    //mean(1) = 100;
    //AdditiveGaussianNoise<2> noise(mean, cov);
    //noise.momentMatrix(2);

    //Covariance<3> cov;
    //cov(0, 0) = 1.00;
    //cov(1, 0) = 6.00;
    //cov(2, 0) = 6.00;
    //cov(0, 1) = 6.00;
    //cov(1, 1) = 4.00;
    //cov(2, 1) = 4.00;
    //cov(0, 2) = 6.00;
    //cov(1, 2) = 4.00;
    //cov(2, 2) = 4.00;

    {
    Eigen::Vector<bry_float_t, 3> mean;
    mean(0) = 10;
    mean(1) = 100;
    mean(2) = 300;
    AdditiveGaussianNoise<3> noise(mean, Covariance<3>::Zero());
    noise.momentMatrix(3);
    }

    //{
    //Eigen::Vector<bry_float_t, 5> mean;
    //mean(0) = 10;
    //mean(1) = 100;
    //mean(2) = 300;
    //mean(3) = 300;
    //mean(4) = 500;
    //AdditiveGaussianNoise<5> noise(mean, Covariance<5>::Zero());
    //noise.momentMatrix(3);
    //}

    //{
    //Eigen::Vector<bry_float_t, 2> mean;
    //mean(0) = 10;
    //mean(1) = 100;
    //AdditiveGaussianNoise<2> noise(mean, Covariance<2>::Zero());
    //noise.momentMatrix(2);
    //}
    return 0;
}
 
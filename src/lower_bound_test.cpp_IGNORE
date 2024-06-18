#include "berry/Polynomial.h"
#include "berry/MultiIndex.h"
#include "berry/Operations.h"
#include "berry/BernsteinTransform.h"

#include "ArgParser.h"
#include "HyperRectangle.h"

#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

using namespace BRY;

int main(int argc, char** argv) {
	ArgParser parser(argc, argv);
	auto deg_inc = parser.parse<bry_int_t>("deg-inc", 'i', 0l, "Degree increase");

    constexpr std::size_t DIM = 2;
    bry_int_t deg = 4;

    Eigen::Tensor<bry_float_t, DIM> tensor(deg + 1, deg + 1);
    tensor.setValues({
        {0.60631369908317678252, 13.38899485314757598076, 18.86898463745490772681, 6.24743038860395927259, 1.34529500295428050549},
        {-2.42525479633270713009, 10.81045211233518443805, 4.97207263775089813862, -13.14568460745550737556, -5.01983046257055320893}, 
        {5.51216124253491912555, -6.43344845860006842031, -36.64165664271996547541, -36.17435354323487217698, -10.50107490507991769846}, 
        {13.15478702354328177648, 12.20662695965430089018, -1.69662035716220316317, -2.04142003212581357730, -0.54996795586680491397}, 
        {2.28913616018212096037, 6.62611192786266656185, 9.88752619335051008420, 6.90385232680833382801, 1.80956659588261459248}});
    tensor = tensor.shuffle(Eigen::array<int, 2>{1, 0});
    //for (bry_int_t i = 0; i < deg; ++i) {
    //    for (bry_int_t j = 0; j < deg; ++j) {
    //        bry_float_t temp = tensor(i, j);
    //        tensor(i, j) = tensor(j, i);
    //        tensor(j, i) = temp;
    //    }
    //}

    Polynomial<DIM> p(tensor);
    //DEBUG("p: " << p);
    Eigen::Map<const Eigen::VectorXd> p_vec(p.tensor().data(), p.nMonomials());
    DEBUG("p vec: " << p_vec.transpose());

    //p.coeff(0, 0) = -.1;
    //p.coeff(1, 0) = -.7;
    //p.coeff(2, 0) = 1.0;
    //p.coeff(0, 1) = -0.4;
    //p.coeff(1, 1) = 0.5;
    //p.coeff(2, 1) = -0.1;
    //p.coeff(0, 2) = 0.5;
    //p.coeff(1, 2) = -0.2;
    //p.coeff(2, 2) = 0.1;


    HyperRectangle<DIM> workspace;
    Eigen::Vector<bry_float_t, DIM> boundary_width{0.2, 0.2};
    workspace.lower_bounds = Eigen::Vector<bry_float_t, DIM>(-1.0, -0.5) - boundary_width;
    workspace.upper_bounds = Eigen::Vector<bry_float_t, DIM>(0.5, 0.5) + boundary_width;

    //workspace.lower_bounds = Eigen::Vector<bry_float_t, DIM>(0.0, 0.0);
    //workspace.upper_bounds = Eigen::Vector<bry_float_t, DIM>(0.9, 0.9);


    //Matrix b_to_p = BernsteinBasisTransform<DIM>::bernToPwrMatrix(deg);
    Matrix p_to_b = BernsteinBasisTransform<DIM>::pwrToBernMatrix(deg, deg_inc.get());

    //DEBUG("ws tf: " << workspace.transformationMatrix(deg));
    Matrix ws_tf = workspace.transformationMatrix(deg);
    Matrix t_mat = p_to_b * ws_tf;
    //DEBUG("tmat: " << t_mat);
    //auto p_bern_tf = transform(p, t_mat);
    Eigen::VectorXd pb_vec_tf = t_mat * p_vec;
    INFO("Lower bound: " << (pb_vec_tf.minCoeff()));

    //DEBUG("tf rows: " << ws_tf.rows() << ", cols: " << ws_tf.cols());
    //NEW_LINE;
    //DEBUG("ws tf: \n" << ws_tf);
    //NEW_LINE;
    //DEBUG("p_vec: " << p_vec.transpose());
    Eigen::VectorXd p_vec_tf = ws_tf * p_vec;
    DEBUG("p_vec tfd: " << p_vec_tf.transpose());
    Polynomial<DIM> p_just_tf(p_vec_tf);

    auto printPCoeffs = [](const Polynomial<DIM>& p) {
        for (std::size_t i = 0; i < p.tensor().size(); ++i) {
            if (i < p.tensor().size() - 1) {
                std::cout << std::fixed << std::setprecision(20) << p.tensor()(i) << ", ";
            } else {
                std::cout << std::fixed << std::setprecision(20) << p.tensor()(i) << std::endl;
            }
        }
    };

    DEBUG("p: ");
    printPCoeffs(p);
    DEBUG("p just tf: ");
    printPCoeffs(p_just_tf);
    return 0;
}

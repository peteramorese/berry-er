#include "Dynamics.h"
#include "Noise.h"
#include "HyperRectangle.h"
#include "Synthesis.h"

#include <iostream>
#include <iomanip>

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

    constexpr std::size_t DIM = 1;
    std::shared_ptr<PolynomialDynamics<DIM>> dynamics_ptr = std::make_shared<PolynomialDynamics<1>>(1);    
    PolynomialDynamics<DIM>& dynamics = *dynamics_ptr;
    dynamics[0].coeff(0) = 0.0;
    dynamics[0].coeff(1) = 0.90;
    //dynamics[0].coeff(2) = 2.2;

    Covariance<DIM> cov;
    cov(0) = 0.01;
    std::shared_ptr<Additive2ndMomentNoise<DIM>> noise_ptr = std::make_shared<Additive2ndMomentNoise<DIM>>(cov);

    // Init set
    std::vector<HyperRectangle<DIM>> init_sets(1);
    init_sets[0].lower_bounds(0) = 0.4;
    init_sets[0].upper_bounds(0) = 0.6;

    // Unsafe set
    std::vector<HyperRectangle<DIM>> unsafe_sets(2);
    unsafe_sets[0].lower_bounds(0) = 0.8;
    unsafe_sets[0].upper_bounds(0) = 1.0;
    unsafe_sets[1].lower_bounds(0) = 0.0;
    unsafe_sets[1].upper_bounds(0) = 0.2;

    // Init set
    std::vector<HyperRectangle<DIM>> safe_sets(1);
    safe_sets[0].lower_bounds(0) = 0.2;
    safe_sets[0].upper_bounds(0) = 0.8;
    //safe_sets[1].lower_bounds(0) = 0.7;
    //safe_sets[1].upper_bounds(0) = 1.0;

    bry_deg_t deg = 10;
    PolyDynamicsSynthesizer synthesizer(dynamics_ptr, noise_ptr, deg);

    synthesizer.setInitialSets(std::move(init_sets));
    synthesizer.setUnsafeSets(std::move(unsafe_sets));
    synthesizer.setSafeSets(std::move(safe_sets));

    synthesizer.initialize();
    auto result = synthesizer.synthesize(5);
    
    INFO("Probability of safety: " << result.p_safe);
    INFO("Eta = " << result.eta << ", Gamma = " << result.gamma);

    Eigen::MatrixXd b_to_p = BernsteinBasisTransform<DIM>::bernToPwrMatrix(deg);
    auto certificate_power = transform(*result.certificate, b_to_p);

    INFO("Barrier: " << certificate_power);
    INFO("Coefficients:");
    for (std::size_t i = 0; i < certificate_power.tensor().size(); ++i) {
        if (i < certificate_power.tensor().size() - 1) {
            std::cout << std::fixed << std::setprecision(20) << certificate_power.tensor()(i) << ", ";
        } else {
            std::cout << std::fixed << std::setprecision(20) << certificate_power.tensor()(i) << std::endl;
        }
    }


    Eigen::MatrixXd phi_inv = BernsteinBasisTransform<DIM>::bernToPwrMatrix(2);
    DEBUG("phi inv: \n" << phi_inv);
    



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


    return 0;
}
 
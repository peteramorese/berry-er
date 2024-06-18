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

    constexpr std::size_t DIM = 2;
    std::shared_ptr<PolynomialDynamics<DIM>> dynamics_ptr = std::make_shared<PolynomialDynamics<DIM>>(1, 1);    
    PolynomialDynamics<DIM>& dynamics = *dynamics_ptr;
    dynamics[0].coeff(0, 0) = 0.01;
    dynamics[0].coeff(1, 0) = 0.90;
    dynamics[0].coeff(0, 1) = 0.90;
    dynamics[0].coeff(1, 1) = 0.95;
    dynamics[1].coeff(0, 0) = 0.01;
    dynamics[1].coeff(1, 0) = 0.90;
    dynamics[1].coeff(0, 1) = 0.90;
    dynamics[1].coeff(1, 1) = 0.95;

    Covariance<DIM> cov;
    cov(0, 0) = 0.01;
    cov(1, 0) = 0.00;
    cov(0, 1) = 0.00;
    cov(1, 1) = 0.01;
    std::shared_ptr<Additive2ndMomentNoise<DIM>> noise_ptr = std::make_shared<Additive2ndMomentNoise<DIM>>(cov);

    std::shared_ptr<SynthesisProblem<DIM>> prob(new SynthesisProblem<DIM>());

    HyperRectangle<DIM> workspace;
    workspace.lower_bounds = Eigen::Vector<bry_float_t, DIM>(-0.1, -0.1);
    workspace.upper_bounds = Eigen::Vector<bry_float_t, DIM>(1.0, 1.0);
    prob->setWorkspace(workspace);

    // Init set
    HyperRectangle<DIM> init_set;
    init_set.lower_bounds(0) = 0.3;
    init_set.upper_bounds(0) = 0.5;
    init_set.lower_bounds(1) = 0.2;
    init_set.upper_bounds(1) = 0.3;
    prob->init_sets.push_back(std::move(init_set));

    // Unsafe set
    HyperRectangle<DIM> unsafe_set_1;
    unsafe_set_1.lower_bounds(0) = 0.2;
    unsafe_set_1.upper_bounds(0) = 0.3;
    unsafe_set_1.lower_bounds(1) = 0.6;
    unsafe_set_1.upper_bounds(1) = 0.7;
    prob->init_sets.push_back(std::move(unsafe_set_1));
    HyperRectangle<DIM> unsafe_set_2;
    unsafe_set_2.lower_bounds(0) = 0.7;
    unsafe_set_2.upper_bounds(0) = 0.8;
    unsafe_set_2.lower_bounds(1) = 0.2;
    unsafe_set_2.upper_bounds(1) = 0.3;
    prob->init_sets.push_back(std::move(unsafe_set_2));

    // Safe set
    prob->safe_sets.resize(7);
    std::vector<HyperRectangle<DIM>>& safe_sets = prob->safe_sets;
    safe_sets[0].lower_bounds(0) = 0.0;
    safe_sets[0].upper_bounds(0) = 0.7;
    safe_sets[0].lower_bounds(1) = 0.0;
    safe_sets[0].upper_bounds(1) = 0.6;

    safe_sets[1].lower_bounds(0) = 0.0;
    safe_sets[1].upper_bounds(0) = 0.2;
    safe_sets[1].lower_bounds(1) = 0.6;
    safe_sets[1].upper_bounds(1) = 0.7;

    safe_sets[2].lower_bounds(0) = 0.0;
    safe_sets[2].upper_bounds(0) = 0.3;
    safe_sets[2].lower_bounds(1) = 0.7;
    safe_sets[2].upper_bounds(1) = 1.0;

    safe_sets[3].lower_bounds(0) = 0.3;
    safe_sets[3].upper_bounds(0) = 0.7;
    safe_sets[3].lower_bounds(1) = 0.6;
    safe_sets[3].upper_bounds(1) = 1.0;

    safe_sets[4].lower_bounds(0) = 0.7;
    safe_sets[4].upper_bounds(0) = 1.0;
    safe_sets[4].lower_bounds(1) = 0.3;
    safe_sets[4].upper_bounds(1) = 1.0;

    safe_sets[5].lower_bounds(0) = 0.8;
    safe_sets[5].upper_bounds(0) = 1.0;
    safe_sets[5].lower_bounds(1) = 0.0;
    safe_sets[5].upper_bounds(1) = 0.3;

    safe_sets[6].lower_bounds(0) = 0.7;
    safe_sets[6].upper_bounds(0) = 0.8;
    safe_sets[6].lower_bounds(1) = 0.0;
    safe_sets[6].upper_bounds(1) = 0.2;

    bry_int_t deg = 10;
    PolyDynamicsSynthesizer synthesizer(dynamics_ptr, noise_ptr, deg);

    prob->time_horizon = 5;

    synthesizer.setProblem(prob);

    synthesizer.setConstraints();
    auto result = synthesizer.synthesize();
    
    INFO("Probability of safety: " << result.p_safe);
    INFO("Eta = " << result.eta << ", Gamma = " << result.gamma);

    //PolynomialDynamics<2> dynamics(1, 1);    
    //dynamics[0].coeff(0, 0) = 5;
    //dynamics[0].coeff(1, 0) = 3;
    //dynamics[0].coeff(1, 1) = 2;

    //dynamics[1].coeff(0, 0) = 6;
    //dynamics[1].coeff(1, 0) = 4;
    //dynamics[1].coeff(1, 1) = 7;

    //DEBUG("Dynamics[0]: " << dynamics[0]);
    //DEBUG("Dynamics[1]: " << dynamics[1]);

    //Matrix dyn_mat = dynamics.dynamicsPowerMatrix(2);
    //Matrix nice_dyn_mat = dyn_mat.unaryExpr([](double x) {
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

    //Matrix add_noise_mat = noise.getAdditiveNoiseMatrix(2);
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
 
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
    DEBUG("double: " << sizeof(double));
    DEBUG("long double: " << sizeof(long double));

    constexpr std::size_t DIM = 1;
    std::shared_ptr<PolynomialDynamics<DIM>> dynamics_ptr = std::make_shared<PolynomialDynamics<1>>(1);    
    PolynomialDynamics<DIM>& dynamics = *dynamics_ptr;
    dynamics[0].coeff(0) = 0.0;
    dynamics[0].coeff(1) = 0.90;

    Covariance<DIM> cov;
    cov(0) = 0.01;
    std::shared_ptr<Additive2ndMomentNoise<DIM>> noise_ptr = std::make_shared<Additive2ndMomentNoise<DIM>>(cov);

    std::shared_ptr<SynthesisProblem<DIM>> prob(new SynthesisProblem<DIM>());

    // Init set
    HyperRectangle<DIM> init_set;
    init_set.lower_bounds(0) = 0.4;
    init_set.upper_bounds(0) = 0.6;
    prob->init_sets.push_back(std::move(init_set));

    // Unsafe sets
    HyperRectangle<DIM> unsafe_set_1;
    unsafe_set_1.lower_bounds(0) = 0.8;
    unsafe_set_1.upper_bounds(0) = 1.0;
    prob->unsafe_sets.push_back(std::move(unsafe_set_1));
    HyperRectangle<DIM> unsafe_set_2;
    unsafe_set_2.lower_bounds(0) = 0.0;
    unsafe_set_2.upper_bounds(0) = 0.2;
    prob->unsafe_sets.push_back(std::move(unsafe_set_2));

    // Safe set
    HyperRectangle<DIM> safe_set;
    safe_set.lower_bounds(0) = 0.2;
    safe_set.upper_bounds(0) = 0.8;
    prob->safe_sets.push_back(std::move(safe_set));

    bry_int_t deg = 10;
    PolyDynamicsSynthesizer synthesizer(dynamics_ptr, noise_ptr, deg);

    prob->time_horizon = 5;
    synthesizer.setProblem(prob);

    synthesizer.initialize();
    auto result = synthesizer.synthesize();
    
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

    return 0;
}
 
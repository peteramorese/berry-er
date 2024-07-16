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

    constexpr std::size_t DIM = 1;

    std::shared_ptr<PolyDynamicsProblem<DIM>> prob(new PolyDynamicsProblem<DIM>());

    prob->dynamics.reset(new PolynomialDynamics<DIM>(2));
    PolynomialDynamics<DIM>& dynamics = *prob->dynamics;
    dynamics[0].coeff(0) = 0.0;
    dynamics[0].coeff(1) = 0.9;
    dynamics[0].coeff(2) = 0.0;
    DEBUG("Dynamics f(x) = " << dynamics[0]);

    Covariance<DIM> cov;
    cov(0, 0) = 0.1;
    prob->noise.reset(new AdditiveGaussianNoise<DIM>(cov));

    // Init set
    HyperRectangle<DIM> init_set;
    init_set.lower_bounds(0) = 0.8;
    init_set.upper_bounds(0) = 0.9;
    prob->init_sets.push_back(init_set);
    
    // Unsafe set
    HyperRectangle<DIM> unsafe_set;
    unsafe_set.lower_bounds(0) = 0.0;
    unsafe_set.upper_bounds(0) = 0.2;
    prob->unsafe_sets.push_back(unsafe_set);

    // Safe set

    HyperRectangle<DIM> safe_set;
    safe_set.lower_bounds(0) = 0.2;
    safe_set.upper_bounds(0) = 1.0;
    prob->safe_sets.push_back(safe_set);

    prob->time_horizon = 2;
    prob->barrier_deg = 3;
    prob->degree_increase = 0;

    prob->subdivide(10);

    Matrix expec_Gamma = prob->noise->additiveNoiseMatrix(prob->barrier_deg);
    DEBUG("E[Gamma] =  \n" << expec_Gamma);
    Matrix F_expec_Gamma = prob->dynamics->dynamicsPowerMatrix(prob->barrier_deg) * expec_Gamma;
    bry_int_t p = prob->dynamics->composedDegree(prob->barrier_deg);
    Matrix Phi_p = BernsteinBasisTransform<DIM>::pwrToBernMatrix(p, 0);
    DEBUG("F expec[Gamma] =  \n" << F_expec_Gamma);

    ConstraintMatrices<DIM> constraints = prob->getConstraintMatrices();

    //INFO("Exporting constraint matrices...");
    writeMatrixToFile(constraints.A, "A.txt");
    writeMatrixToFile(constraints.b, "b.txt");
    writeMatrixToFile(Phi_p, "Phi_p.txt");
    writeMatrixToFile(safe_set.transformationMatrix(p), "T.txt");
    //INFO("Done!");

    //if (solve) {
        INFO("Solving...");
    //    
        auto result = synthesize(constraints, prob->time_horizon, "gurobi");
        INFO("Done!");
    //    NEW_LINE;
        INFO("Probability of safety: " << result.p_safe);

        printf("Eta = %.32f\n", result.eta);
        printf("Gamma = %.32f\n", result.gamma);
        DEBUG("b values: " << result.b_values.transpose());
    DEBUG("Diff coeffs: " << ((F_expec_Gamma - Matrix::Identity(F_expec_Gamma.rows(), F_expec_Gamma.cols()))*result.b_values).transpose());
    //DEBUG("Lower bounds: " << 
    //    (
    //    -Phi_p*safe_set.transformationMatrix(p)*(F_expec_Gamma - Matrix::Identity(F_expec_Gamma.rows(), F_expec_Gamma.cols()))*result.b_values + result.gamma*Vector::Ones(Phi_p.rows())
    //    ).transpose());

    
    //    //INFO("Eta = " << result.eta << ", Gamma = " << result.gamma);
    //    INFO("Computation time: " << result.comp_time << "s");

    //    if (result.diagDeg())
    //        result.fromDiagonalDegree();

    //    writeMatrixToFile(result.b_values, "certificate_coeffs.txt");
    //}


    return 0;
}
 
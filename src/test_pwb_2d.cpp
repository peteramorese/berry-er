#include "Dynamics.h"
#include "Noise.h"
#include "HyperRectangle.h"
#include "Synthesis.h"
#include "ArgParser.h"

#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

using namespace BRY;

int main(int argc, char** argv) {

	ArgParser parser(argc, argv);
	bool cp_coeffs = parser.parse<void>("cp-coeffs", 'c', "Print the copy-pastable coefficients").has();
	auto solver_id = parser.parse<std::string>("s-id", 's', "SCIP", "Solver ID");
	auto barrier_deg = parser.parse<bry_int_t>("deg", 'd', 4l, "Barrier degree");
	auto deg_increase = parser.parse<bry_int_t>("deg-inc", 'i', 0l, "Barrier degree increase");
	auto time_steps = parser.parse<uint64_t>("ts", 't', 5, "Number of time steps");
    parser.enableHelp();

    constexpr std::size_t DIM = 2;
    std::shared_ptr<PolynomialDynamics<DIM>> dynamics_ptr = std::make_shared<PolynomialDynamics<DIM>>(1, 1);    
    PolynomialDynamics<DIM>& dynamics = *dynamics_ptr;
    dynamics[0].coeff(0, 0) = 0.0;
    dynamics[0].coeff(1, 0) = 0.5;
    dynamics[0].coeff(0, 1) = 0.0;
    dynamics[0].coeff(1, 1) = 0.0;
    dynamics[1].coeff(0, 0) = 0.0;
    dynamics[1].coeff(1, 0) = 0.0;
    dynamics[1].coeff(0, 1) = 0.5;
    dynamics[1].coeff(1, 1) = 0.0;

    Covariance<DIM> cov;
    cov(0, 0) = 0.01;
    cov(1, 0) = 0.00;
    cov(0, 1) = 0.00;
    cov(1, 1) = 0.01;
    std::shared_ptr<Additive2ndMomentNoise<DIM>> noise_ptr = std::make_shared<Additive2ndMomentNoise<DIM>>(cov);

    std::shared_ptr<SynthesisProblem<DIM>> prob(new SynthesisProblem<DIM>());

    Eigen::Vector<bry_float_t, DIM> boundary_width{0.2, 0.2};

    prob->workspace.lower_bounds = Eigen::Vector<bry_float_t, DIM>(-1.0, -0.5) - boundary_width;
    prob->workspace.upper_bounds = Eigen::Vector<bry_float_t, DIM>(0.5, 0.5) + boundary_width;

    auto printSetBounds = [](const HyperRectangle<DIM>& set) {
        DEBUG("Set bounds: [" 
            << set.lower_bounds(0) << ", " 
            << set.upper_bounds(0) << ", "
            << set.lower_bounds(1) << ", " 
            << set.upper_bounds(1) << "]");
    };

    // Init set
    HyperRectangle<DIM> init_set;
    init_set.lower_bounds(0) = -0.8;
    init_set.upper_bounds(0) = -0.6;
    init_set.lower_bounds(1) = -0.2;
    init_set.upper_bounds(1) = 0.0;
    prob->init_sets.push_back(init_set);
    DEBUG("Init set:");
    printSetBounds(init_set);
    NEW_LINE;

    DEBUG("Unsafe sets:");
    // Unsafe set
    HyperRectangle<DIM> boundary_left;
    // Boundary left
    boundary_left.lower_bounds(0) = -1.0 - boundary_width(0);
    boundary_left.upper_bounds(0) = -1.0;
    boundary_left.lower_bounds(1) = -0.5 - boundary_width(1);
    boundary_left.upper_bounds(1) = 0.5 + boundary_width(1);
    prob->unsafe_sets.push_back(boundary_left);
    printSetBounds(boundary_left);
    // Boundary right
    HyperRectangle<DIM> boundary_right;
    boundary_right.lower_bounds(0) = 0.5;
    boundary_right.upper_bounds(0) = 0.5 + boundary_width(0);
    boundary_right.lower_bounds(1) = -0.5 - boundary_width(1);
    boundary_right.upper_bounds(1) = 0.5 + boundary_width(1);
    prob->unsafe_sets.push_back(boundary_right);
    printSetBounds(boundary_right);
    // Boundary top
    HyperRectangle<DIM> boundary_top;
    boundary_top.lower_bounds(0) = -1.0;
    boundary_top.upper_bounds(0) = 0.5;
    boundary_top.lower_bounds(1) = 0.5;
    boundary_top.upper_bounds(1) = 0.5 + boundary_width(1);
    prob->unsafe_sets.push_back(boundary_top);
    printSetBounds(boundary_top);
    // Boundary bottom
    HyperRectangle<DIM> boundary_bottom;
    boundary_bottom.lower_bounds(0) = -1.0;
    boundary_bottom.upper_bounds(0) = 0.5;
    boundary_bottom.lower_bounds(1) = -0.5 - boundary_width(1);
    boundary_bottom.upper_bounds(1) = -0.5;
    prob->unsafe_sets.push_back(boundary_bottom);
    printSetBounds(boundary_bottom);
    // Non convex unsafe regions
    HyperRectangle<DIM> upper_region;
    upper_region.lower_bounds(0) = -0.57;
    upper_region.upper_bounds(0) = -0.53;
    upper_region.lower_bounds(1) = -0.17;
    upper_region.upper_bounds(1) = -0.13;
    prob->unsafe_sets.push_back(upper_region);
    printSetBounds(upper_region);
    HyperRectangle<DIM> lower_region;
    lower_region.lower_bounds(0) = -0.57;
    lower_region.upper_bounds(0) = -0.53;
    lower_region.lower_bounds(1) = 0.28;
    lower_region.upper_bounds(1) = 0.32;
    prob->unsafe_sets.push_back(lower_region);
    printSetBounds(lower_region);

    // Safe set
    HyperRectangle<DIM> safe_set;
    safe_set.lower_bounds(0) = -1.0;
    safe_set.upper_bounds(0) = 0.5;
    safe_set.lower_bounds(1) = -0.5;
    safe_set.upper_bounds(1) = 0.5;
    prob->safe_sets.push_back(safe_set);

    //prob->safe_sets.resize(5);
    //std::vector<HyperRectangle<DIM>>& safe_sets = prob->safe_sets;
    //safe_sets[0].lower_bounds(0) = -1.0;
    //safe_sets[0].upper_bounds(0) = -0.57;
    //safe_sets[0].lower_bounds(1) = -0.5;
    //safe_sets[0].upper_bounds(1) = 0.5;

    //safe_sets[1].lower_bounds(0) = -0.57;
    //safe_sets[1].upper_bounds(0) = -0.53;
    //safe_sets[1].lower_bounds(1) = -0.5;
    //safe_sets[1].upper_bounds(1) = -0.17;

    //safe_sets[2].lower_bounds(0) = -0.57;
    //safe_sets[2].upper_bounds(0) = -0.53;
    //safe_sets[2].lower_bounds(1) = -0.13;
    //safe_sets[2].upper_bounds(1) = 0.28;

    //safe_sets[3].lower_bounds(0) = -0.57;
    //safe_sets[3].upper_bounds(0) = -0.53;
    //safe_sets[3].lower_bounds(1) = 0.32;
    //safe_sets[3].upper_bounds(1) = 0.5;

    //safe_sets[4].lower_bounds(0) = -0.53;
    //safe_sets[4].upper_bounds(0) = 0.5;
    //safe_sets[4].lower_bounds(1) = -0.5;
    //safe_sets[4].upper_bounds(1) = 0.5;


    //bry_int_t deg = 8;
    PolyDynamicsSynthesizer synthesizer(dynamics_ptr, noise_ptr, barrier_deg.get(), solver_id.get());

    prob->time_horizon = time_steps.get();
    synthesizer.setProblem(prob);

    synthesizer.initialize(deg_increase.get());
    auto result = synthesizer.synthesize();
    
    INFO("Probability of safety: " << result.p_safe);
    INFO("Eta = " << result.eta << ", Gamma = " << result.gamma);

    Matrix b_to_p = BernsteinBasisTransform<DIM>::bernToPwrMatrix(barrier_deg.get());
    auto certificate_power = transform(*result.certificate, b_to_p);

    INFO("Barrier: " << certificate_power);

    if (cp_coeffs) {
        NEW_LINE;
        for (std::size_t i = 0; i < certificate_power.tensor().size(); ++i) {
            if (i < certificate_power.tensor().size() - 1) {
                std::cout << std::fixed << std::setprecision(20) << certificate_power.tensor()(i) << ", ";
            } else {
                std::cout << std::fixed << std::setprecision(20) << certificate_power.tensor()(i) << std::endl;
            }
        }
        NEW_LINE;
    }

    return 0;
}
 
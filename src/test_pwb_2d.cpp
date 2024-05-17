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
	auto barrier_deg = parser.parse<bry_deg_t>("deg", 'd', 4l, "Barrier degree");
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

    Eigen::Vector<bry_float_t, DIM> boundary_width{0.2, 0.2};

    HyperRectangle<DIM> workspace;
    workspace.lower_bounds = Eigen::Vector<bry_float_t, DIM>(-1.0, -0.5) - boundary_width;
    workspace.upper_bounds = Eigen::Vector<bry_float_t, DIM>(0.5, 0.5) + boundary_width;

    auto printSetBounds = [](const HyperRectangle<DIM>& set) {
        DEBUG("Set bounds: [" 
            << set.lower_bounds(0) << ", " 
            << set.upper_bounds(0) << ", "
            << set.lower_bounds(1) << ", " 
            << set.upper_bounds(1) << "]");
    };

    // Init set
    std::vector<HyperRectangle<DIM>> init_sets(1);
    init_sets[0].lower_bounds(0) = -0.8;
    init_sets[0].upper_bounds(0) = -0.6;
    init_sets[0].lower_bounds(1) = -0.2;
    init_sets[0].upper_bounds(1) = 0.0;
    DEBUG("Init set:");
    printSetBounds(init_sets[0]);
    NEW_LINE;

    DEBUG("Unsafe sets:");
    // Unsafe set
    std::vector<HyperRectangle<DIM>> unsafe_sets(4);
    // Boundary left
    unsafe_sets[0].lower_bounds(0) = -1.0 - boundary_width(0);
    unsafe_sets[0].upper_bounds(0) = -1.0;
    unsafe_sets[0].lower_bounds(1) = -0.5 - boundary_width(1);
    unsafe_sets[0].upper_bounds(1) = 0.5 + boundary_width(1);
    printSetBounds(unsafe_sets[0]);
    // Boundary right
    unsafe_sets[1].lower_bounds(0) = 0.5;
    unsafe_sets[1].upper_bounds(0) = 0.5 + boundary_width(0);
    unsafe_sets[1].lower_bounds(1) = -0.5 - boundary_width(1);
    unsafe_sets[1].upper_bounds(1) = 0.5 + boundary_width(1);
    printSetBounds(unsafe_sets[1]);
    // Boundary top
    unsafe_sets[2].lower_bounds(0) = -1.0;
    unsafe_sets[2].upper_bounds(0) = 0.5;
    unsafe_sets[2].lower_bounds(1) = 0.5;
    unsafe_sets[2].upper_bounds(1) = 0.5 + boundary_width(1);
    printSetBounds(unsafe_sets[2]);
    // Boundary bottom
    unsafe_sets[3].lower_bounds(0) = -1.0;
    unsafe_sets[3].upper_bounds(0) = 0.5;
    unsafe_sets[3].lower_bounds(1) = -0.5 - boundary_width(1);
    unsafe_sets[3].upper_bounds(1) = -0.5;
    printSetBounds(unsafe_sets[3]);
    // Non convex unsafe regions
    //unsafe_sets[4].lower_bounds(0) = -0.57;
    //unsafe_sets[4].upper_bounds(0) = -0.53;
    //unsafe_sets[4].lower_bounds(1) = -0.17;
    //unsafe_sets[4].upper_bounds(1) = -0.13;
    //unsafe_sets[5].lower_bounds(0) = -0.57;
    //unsafe_sets[5].upper_bounds(0) = -0.53;
    //unsafe_sets[5].lower_bounds(1) = 0.28;
    //unsafe_sets[5].upper_bounds(1) = 0.32;

    // Safe set
    std::vector<HyperRectangle<DIM>> safe_sets(5);
    safe_sets[0].lower_bounds(0) = -1.0;
    safe_sets[0].upper_bounds(0) = -0.57;
    safe_sets[0].lower_bounds(1) = -0.5;
    safe_sets[0].upper_bounds(1) = 0.5;

    safe_sets[1].lower_bounds(0) = -0.57;
    safe_sets[1].upper_bounds(0) = -0.53;
    safe_sets[1].lower_bounds(1) = -0.5;
    safe_sets[1].upper_bounds(1) = -0.17;

    safe_sets[2].lower_bounds(0) = -0.57;
    safe_sets[2].upper_bounds(0) = -0.53;
    safe_sets[2].lower_bounds(1) = -0.13;
    safe_sets[2].upper_bounds(1) = 0.28;

    safe_sets[3].lower_bounds(0) = -0.57;
    safe_sets[3].upper_bounds(0) = -0.53;
    safe_sets[3].lower_bounds(1) = 0.32;
    safe_sets[3].upper_bounds(1) = 0.5;

    safe_sets[4].lower_bounds(0) = -0.53;
    safe_sets[4].upper_bounds(0) = 0.5;
    safe_sets[4].lower_bounds(1) = -0.5;
    safe_sets[4].upper_bounds(1) = 0.5;


    //bry_deg_t deg = 8;
    PolyDynamicsSynthesizer synthesizer(dynamics_ptr, noise_ptr, barrier_deg.get(), solver_id.get());

    synthesizer.setWorkspace(workspace);
    synthesizer.setInitialSets(std::move(init_sets));
    synthesizer.setUnsafeSets(std::move(unsafe_sets));
    synthesizer.setSafeSets(std::move(safe_sets));

    synthesizer.initialize();
    auto result = synthesizer.synthesize(time_steps.get());
    
    INFO("Probability of safety: " << result.p_safe);
    INFO("Eta = " << result.eta << ", Gamma = " << result.gamma);

    Eigen::MatrixXd b_to_p = BernsteinBasisTransform<DIM>::bernToPwrMatrix(barrier_deg.get());
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
 
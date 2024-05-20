#include "Dynamics.h"
#include "Noise.h"
#include "HyperRectangle.h"
#include "Synthesis.h"
#include "ArgParser.h"

#include <iostream>
#include <fstream>
#include <iomanip>

#include <Eigen/Dense>

using namespace BRY;

void writeMatrixToFile(const Eigen::MatrixXd& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        // Write matrix dimensions
        file << matrix.rows() << " " << matrix.cols() << std::endl;

        // Write matrix data
        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = 0; j < matrix.cols(); ++j) {
                file << std::fixed << std::setprecision(20) << matrix(i, j);
                if (j < matrix.cols() - 1) file << " ";
            }
            file << std::endl;
        }

        file.close();
    } else {
        std::cerr << "Unable to open file";
    }
}

//#define SIMPLE_PROBLEM

int main(int argc, char** argv) {

	ArgParser parser(argc, argv);
    bool verbose = parser.parse<void>('v', "Verbose").has();
    bool solve = parser.parse<void>('r', "Solve the synthesis problem").has();
	auto solver_id = parser.parse<std::string>("s-id", 's', "SCIP", "Solver ID");
	auto barrier_deg = parser.parse<bry_deg_t>("deg", 'd', 4l, "Barrier degree");
	auto time_steps = parser.parse<uint64_t>("ts", 't', 5, "Number of time steps");
    parser.enableHelp();

#ifdef SIMPLE_PROBLEM
    constexpr std::size_t DIM = 1;
    std::shared_ptr<PolynomialDynamics<DIM>> dynamics_ptr = std::make_shared<PolynomialDynamics<1>>(1);    
    PolynomialDynamics<DIM>& dynamics = *dynamics_ptr;
    dynamics[0].coeff(0) = 0.0;
    dynamics[0].coeff(1) = 0.90;
    //dynamics[0].coeff(2) = 2.2;

    Covariance<DIM> cov;
    cov(0) = 0.01;
    std::shared_ptr<Additive2ndMomentNoise<DIM>> noise_ptr = std::make_shared<Additive2ndMomentNoise<DIM>>(cov);

    HyperRectangle<DIM> workspace;

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
#else
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
    std::vector<HyperRectangle<DIM>> safe_sets(1);
    safe_sets[0].lower_bounds(0) = -1.0;
    safe_sets[0].upper_bounds(0) = 0.5;
    safe_sets[0].lower_bounds(1) = -0.5;
    safe_sets[0].upper_bounds(1) = 0.5;

    //std::vector<HyperRectangle<DIM>> safe_sets(5);
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
#endif

    PolyDynamicsSynthesizer synthesizer(dynamics_ptr, noise_ptr, barrier_deg.get(), solver_id.get());

    synthesizer.setWorkspace(workspace);
    synthesizer.setInitialSets(std::move(init_sets));
    synthesizer.setUnsafeSets(std::move(unsafe_sets));
    synthesizer.setSafeSets(std::move(safe_sets));

    Eigen::MatrixXd b_to_p = BernsteinBasisTransform<DIM>::bernToPwrMatrix(barrier_deg.get());

    auto[A, b] = synthesizer.getConstraintMatrices();
    if (verbose) {
        INFO("A: \n" << A);
        NEW_LINE;
        INFO("b: \n" << b.transpose());
    }
    writeMatrixToFile(A, "A.txt");
    writeMatrixToFile(b, "b.txt");
    writeMatrixToFile(b_to_p, "Phi_inv.txt");
    if (solve) {
        synthesizer.initialize();
        auto result = synthesizer.synthesize(time_steps.get());
        
        INFO("Probability of safety: " << result.p_safe);
        INFO("Eta = " << result.eta << ", Gamma = " << result.gamma);

        auto certificate_power = transform(*result.certificate, b_to_p);

        INFO("Barrier: " << certificate_power);
    }



    //INFO("Coefficients:");
    //for (std::size_t i = 0; i < certificate_power.tensor().size(); ++i) {
    //    if (i < certificate_power.tensor().size() - 1) {
    //        std::cout << std::fixed << std::setprecision(20) << certificate_power.tensor()(i) << ", ";
    //    } else {
    //        std::cout << std::fixed << std::setprecision(20) << certificate_power.tensor()(i) << std::endl;
    //    }
    //}


    //Eigen::MatrixXd phi_inv = BernsteinBasisTransform<DIM>::bernToPwrMatrix(2);
    //DEBUG("phi inv: \n" << phi_inv);
    



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
 
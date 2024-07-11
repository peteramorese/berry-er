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

	ArgParser parser(argc, argv);
    bool verbose = parser.parse<void>('v', "Verbose").has();
    bool solve = parser.parse<void>('r', "Solve the synthesis problem").has();
    bool non_convex = parser.parse<void>("non-conv", "Solve the non-convex synthesis problem (default to convex)").has();
    bool diag_deg = parser.parse<void>("diag-deg", "Use the diagonal polynomial degree definition").has();
	auto solver_id = parser.parse<std::string>("s-id", 's', "SCIP", "Solver ID");
	auto barrier_deg = parser.parse<bry_int_t>("deg", 'd', 4l, "Barrier degree");
	auto deg_increase = parser.parse<bry_int_t>("deg-inc", 'i', 0l, "Barrier degree increase");
	auto subd = parser.parse<bry_int_t>("subdiv", "Set subdivision");
	auto time_steps = parser.parse<uint64_t>("ts", 't', 10, "Number of time steps");
    parser.enableHelp();

    constexpr std::size_t DIM = 2;

    std::shared_ptr<PolyDynamicsProblem<DIM>> prob(new PolyDynamicsProblem<DIM>());

    prob->dynamics.reset(new PolynomialDynamics<DIM>(1, 1));
    PolynomialDynamics<DIM>& dynamics = *prob->dynamics;
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
    prob->noise.reset(new AdditiveGaussianNoise<DIM>(cov));

    Eigen::Vector<bry_float_t, DIM> boundary_width{0.2, 0.2};

    auto printSetBounds = [](const HyperRectangle<DIM>& set) {
        DEBUG("Set bounds: [" 
            << set.lower_bounds(0) << ", " 
            << set.upper_bounds(0) << ", "
            << set.lower_bounds(1) << ", " 
            << set.upper_bounds(1) << "]");
    };


    HyperRectangle<DIM> workspace;
    workspace.lower_bounds = Eigen::Vector<bry_float_t, DIM>(-1.0, -0.5) - boundary_width;
    workspace.upper_bounds = Eigen::Vector<bry_float_t, DIM>(0.5, 0.5) + boundary_width;
    prob->setWorkspace(workspace);
    //DEBUG("Workspace set:");
    //printSetBounds(workspace);

    // Init set
    HyperRectangle<DIM> init_set;
    init_set.lower_bounds(0) = -0.8;
    init_set.upper_bounds(0) = -0.6;
    init_set.lower_bounds(1) = 0.0;
    init_set.upper_bounds(1) = 0.2;
    prob->init_sets.push_back(init_set);
    //DEBUG("Init set:");
    //printSetBounds(init_set);
    
    NEW_LINE;

    //DEBUG("Unsafe sets:");
    // Unsafe set
    HyperRectangle<DIM> boundary_left;
    // Boundary left
    boundary_left.lower_bounds(0) = -1.0 - boundary_width(0);
    boundary_left.upper_bounds(0) = -1.0;
    boundary_left.lower_bounds(1) = -0.5 - boundary_width(1);
    boundary_left.upper_bounds(1) = 0.5 + boundary_width(1);
    prob->unsafe_sets.push_back(boundary_left);
    //printSetBounds(boundary_left);
    // Boundary right
    HyperRectangle<DIM> boundary_right;
    boundary_right.lower_bounds(0) = 0.5;
    boundary_right.upper_bounds(0) = 0.5 + boundary_width(0);
    boundary_right.lower_bounds(1) = -0.5 - boundary_width(1);
    boundary_right.upper_bounds(1) = 0.5 + boundary_width(1);
    prob->unsafe_sets.push_back(boundary_right);
    //printSetBounds(boundary_right);
    // Boundary top
    HyperRectangle<DIM> boundary_top;
    boundary_top.lower_bounds(0) = -1.0;
    boundary_top.upper_bounds(0) = 0.5;
    boundary_top.lower_bounds(1) = 0.5;
    boundary_top.upper_bounds(1) = 0.5 + boundary_width(1);
    prob->unsafe_sets.push_back(boundary_top);
    //printSetBounds(boundary_top);
    // Boundary bottom
    HyperRectangle<DIM> boundary_bottom;
    boundary_bottom.lower_bounds(0) = -1.0;
    boundary_bottom.upper_bounds(0) = 0.5;
    boundary_bottom.lower_bounds(1) = -0.5 - boundary_width(1);
    boundary_bottom.upper_bounds(1) = -0.5;
    prob->unsafe_sets.push_back(boundary_bottom);
    //printSetBounds(boundary_bottom);
    if (non_convex) {
        // Non convex unsafe regions
        HyperRectangle<DIM> upper_region;
        upper_region.lower_bounds(0) = -0.57;
        upper_region.upper_bounds(0) = -0.53;
        upper_region.lower_bounds(1) = -0.17;
        upper_region.upper_bounds(1) = -0.13;
        prob->unsafe_sets.push_back(upper_region);
        //printSetBounds(upper_region);
        HyperRectangle<DIM> lower_region;
        lower_region.lower_bounds(0) = -0.57;
        lower_region.upper_bounds(0) = -0.53;
        lower_region.lower_bounds(1) = 0.28;
        lower_region.upper_bounds(1) = 0.32;
        prob->unsafe_sets.push_back(lower_region);
        //printSetBounds(lower_region);
    }

    // Safe set

    if (!non_convex) {
        HyperRectangle<DIM> safe_set;
        safe_set.lower_bounds(0) = -1.0;
        safe_set.upper_bounds(0) = 0.5;
        safe_set.lower_bounds(1) = -0.5;
        safe_set.upper_bounds(1) = 0.5;
        prob->safe_sets.push_back(safe_set);
    } else {
        std::list<HyperRectangle<DIM>>& safe_sets = prob->safe_sets;
        {
            HyperRectangle<DIM> set;
            set.lower_bounds(0) = -1.0;
            set.upper_bounds(0) = -0.57;
            set.lower_bounds(1) = -0.5;
            set.upper_bounds(1) = 0.5;
            safe_sets.push_back(set);
        }

        {
            HyperRectangle<DIM> set;
            set.lower_bounds(0) = -0.57;
            set.upper_bounds(0) = -0.53;
            set.lower_bounds(1) = -0.5;
            set.upper_bounds(1) = -0.17;
            safe_sets.push_back(set);
        }

        {
            HyperRectangle<DIM> set;
            set.lower_bounds(0) = -0.57;
            set.upper_bounds(0) = -0.53;
            set.lower_bounds(1) = -0.13;
            set.upper_bounds(1) = 0.28;
            safe_sets.push_back(set);
        }

        {
            HyperRectangle<DIM> set;
            set.lower_bounds(0) = -0.57;
            set.upper_bounds(0) = -0.53;
            set.lower_bounds(1) = 0.32;
            set.upper_bounds(1) = 0.5;
            safe_sets.push_back(set);
        }

        {
            HyperRectangle<DIM> set;
            set.lower_bounds(0) = -0.53;
            set.upper_bounds(0) = 0.5;
            set.lower_bounds(1) = -0.5;
            set.upper_bounds(1) = 0.5;
            safe_sets.push_back(set);
        }
    }

    prob->time_horizon = time_steps.get();
    prob->barrier_deg = barrier_deg.get();
    prob->degree_increase = deg_increase.get();
    prob->diag_deg = diag_deg;

    if (subd.has()) {
        INFO("Subdividing in " << subd.get());
        prob->subdivide(subd.get());
    }

    //Matrix F_expec_Gamma = prob->dynamics->dynamicsPowerMatrix(prob->barrier_deg) * prob->noise->additiveNoiseMatrix(prob->barrier_deg);
    //DEBUG("F expec[Gamma] =  \n" << F_expec_Gamma);
    //writeMatrixToFile(F_expec_Gamma, "F_expec_Gamma.txt");

    //auto typeToStr = [](ConstraintType type) {
    //    std::string s;
    //    switch (type) {
    //        case ConstraintType::Workspace:
    //            s = "Workspace";
    //            break;
    //        case ConstraintType::Init:
    //            s = "Init";
    //            break;
    //        case ConstraintType::Unsafe:
    //            s = "Unsafe";
    //            break;
    //        case ConstraintType::Safe:
    //            s = "Safe";
    //            break;
    //    }
    //    return s;
    //};

    ConstraintMatrices<DIM> constraints = prob->getConstraintMatrices();

    //int i = 0;
    //for (auto[type, idx] : constraints.constraint_ids) {
    //    DEBUG("Constraint " << i++ << ", ID type: " << typeToStr(type) << ", set idx: " << idx);
    //}

    //bry_int_t idx = 1194;
    //auto set = *prob->lookupSetFromConstraint(constraints.constraint_ids[idx]);
    //DEBUG("Looked up constraint: " << idx << " and got set: ");
    //printSetBounds(set);

    INFO("Exporting constraint matrices...");
    writeMatrixToFile(constraints.A, "A.txt");
    writeMatrixToFile(constraints.b, "b.txt");
    INFO("Done!");

    if (solve) {
        INFO("Solving...");
        
        auto result = synthesize(constraints, prob->time_horizon, solver_id.get());
        INFO("Done!");
        NEW_LINE;
        INFO("Probability of safety: " << result.p_safe);

        printf("Eta = %.32f\n", result.eta);
        printf("Gamma = %.32f\n", result.gamma);
        //INFO("Eta = " << result.eta << ", Gamma = " << result.gamma);
        INFO("Computation time: " << result.comp_time << "s");

        if (result.diagDeg())
            result.fromDiagonalDegree();

        writeMatrixToFile(result.b_values, "certificate_coeffs.txt");
    }


    return 0;
}
 
#include "Dynamics.h"
#include "Noise.h"
#include "HyperRectangle.h"
#include "Synthesis.h"

#include "lemon/ArgParser.h"

#include <iostream>
#include <stdio.h>

#include <Eigen/Dense>

using namespace BRY;

int main(int argc, char** argv) {

	lemon::ArgParser parser(argc, argv);
    lemon::Arg<lemon::ArgT::Check> verbose = parser.addDef<lemon::ArgT::Check>().flag('v').key("Verbose");
    lemon::Arg<lemon::ArgT::Check> non_convex = parser.addDef<lemon::ArgT::Check>().key("non-conv").description("Solve the non-convex synthesis problem (default to convex)");
    lemon::Arg<lemon::ArgT::Check> adaptive = parser.addDef<lemon::ArgT::Check>().flag('a').description("Use the adaptive subdivision algorithm");
    lemon::Arg<lemon::ArgT::Check> export_matrices = parser.addDef<lemon::ArgT::Check>().flag('e').description("Export the matrices to use external solvers");
    lemon::Arg<lemon::ArgT::Value, std::string> filter = parser.addDef<lemon::ArgT::Value, std::string>().flag('f').key("filter").description("Select which filter to use").options({"diagdeg", "oddsum"});
    lemon::Arg<lemon::ArgT::Value, std::string> solver_id = parser.addDef<lemon::ArgT::Value, std::string>().flag('s').key("s-id").description("Solver ID").defaultValue("SCIP");
    lemon::Arg<lemon::ArgT::Value, bry_int_t> barrier_deg = parser.addDef<lemon::ArgT::Value, bry_int_t>().flag('d').key("deg").description("Barrier degree").required();
	lemon::Arg<lemon::ArgT::Value, bry_int_t> deg_increase = parser.addDef<lemon::ArgT::Value, bry_int_t>().flag('i').key("deg-inc").defaultValue(0l).description("Barrier degree increase");
	lemon::Arg<lemon::ArgT::Value, bry_int_t> subd = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("subd").flag('s').defaultValue(0l).description("Barrier subdivision");
	lemon::Arg<lemon::ArgT::Value, bry_int_t> ada_iters = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("ada-iters").defaultValue(1l).description("Number of adaptive subdivision iterations (ONLY FOR ADAPTIVE)");
	lemon::Arg<lemon::ArgT::Value, bry_int_t> ada_max_subdiv = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("ada-max-subdiv").defaultValue(2l).description("Max number of sets divided each iteration (ONLY FOR ADAPTIVE)");
	lemon::Arg<lemon::ArgT::Value, bry_int_t> time_steps = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("ts").flag('t').defaultValue(10l).description("Number of time steps");
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

    prob->time_horizon = time_steps.value();
    prob->barrier_deg = barrier_deg.value();
    prob->degree_increase = deg_increase.value();
    if (filter) {
        if (filter.value() == "diagdeg") {
            prob->filter = std::make_shared<DiagDegFilter<DIM>>(barrier_deg.value());
        } else if (filter.value() == "oddsum") {
            prob->filter = std::make_shared<OddSumFilter<DIM>>();
        }
    }

    if (subd) {
        INFO("Subdividing in " << subd.value());
        prob->subdivide(subd.value());
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


    //int i = 0;
    //for (auto[type, idx] : constraints.constraint_ids) {
    //    DEBUG("Constraint " << i++ << ", ID type: " << typeToStr(type) << ", set idx: " << idx);
    //}

    //bry_int_t idx = 1194;
    //auto set = *prob->lookupSetFromConstraint(constraints.constraint_ids[idx]);
    //DEBUG("Looked up constraint: " << idx << " and got set: ");
    //printSetBounds(set);

    if (export_matrices) {
        ConstraintMatrices<DIM> constraints = prob->getConstraintMatrices();
        INFO("Exporting constraint matrices...");
        writeMatrixToFile(constraints.A, "A.txt");
        writeMatrixToFile(constraints.b, "b.txt");
        INFO("Done!");
    }

    INFO("Solving...");
    SynthesisResult<DIM> result;
    if (adaptive) {
        result = synthesizeAdaptive(*prob, ada_iters.value(), ada_max_subdiv.value(), solver_id.value());
    } else {
        result = synthesize(*prob, solver_id.value());
    }
    INFO("Done!");
    NEW_LINE;
    INFO("Probability of safety: " << result.p_safe);

    printf("Eta = %.32f\n", result.eta);
    printf("Gamma = %.32f\n", result.gamma);
    //INFO("Eta = " << result.eta << ", Gamma = " << result.gamma);
    INFO("Computation time: " << result.comp_time << "s");

    if (result.isFilterApplied()) {
        result.removeFilter();
    }

    writeMatrixToFile(result.b_values, "certificate_coeffs.txt");


    return 0;
}
 
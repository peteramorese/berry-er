#include "Dynamics.h"
#include "Noise.h"
#include "HyperRectangle.h"
#include "Synthesis.h"
#include "AdaptiveProblem.h"
#include "Tools.h"

#include "lemon/ArgParser.h"

#include <iostream>
#include <stdio.h>

#include <Eigen/Dense>

using namespace BRY;

constexpr std::size_t DIM = 2;

int main(int argc, char** argv) {

	lemon::ArgParser parser(argc, argv);
    lemon::Arg<lemon::ArgT::Check> verbose = parser.addDef<lemon::ArgT::Check>().flag('v').key("Verbose");
    lemon::Arg<lemon::ArgT::Check> non_convex = parser.addDef<lemon::ArgT::Check>().key("non-conv").description("Solve the non-convex synthesis problem (default to convex)");
    //lemon::Arg<lemon::ArgT::Check> adaptive = parser.addDef<lemon::ArgT::Check>().flag('a').description("Use the adaptive subdivision algorithm");
    lemon::Arg<lemon::ArgT::Check> export_matrices = parser.addDef<lemon::ArgT::Check>().flag('e').description("Export the matrices to use external solvers");
    lemon::Arg<lemon::ArgT::Value, std::string> filter = parser.addDef<lemon::ArgT::Value, std::string>().flag('f').key("filter").description("Select which filter to use").options({"diagdeg", "oddsum"});
    lemon::Arg<lemon::ArgT::Value, std::string> solver_id = parser.addDef<lemon::ArgT::Value, std::string>().key("solver").description("Solver ID").defaultValue("SCIP");
    lemon::Arg<lemon::ArgT::Value, std::string> dynamics_type = parser.addDef<lemon::ArgT::Value, std::string>().key("dynamics-type").description("Type of dynamics").defaultValue("to_origin").options({"to_origin", "random"});
	lemon::Arg<lemon::ArgT::Value, bry_int_t> dynamics_deg = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("dynamics-deg").defaultValue(1l).description("Degree of dynamics (e.g. 1 is linear, 2 is quadratic, etc.) (ONLY FOR `random` DYNAMICS)");
    lemon::Arg<lemon::ArgT::Value, bry_int_t> barrier_deg = parser.addDef<lemon::ArgT::Value, bry_int_t>().flag('d').key("deg").description("Barrier degree").required();
	lemon::Arg<lemon::ArgT::Value, bry_int_t> deg_increase = parser.addDef<lemon::ArgT::Value, bry_int_t>().flag('i').key("deg-inc").defaultValue(0l).description("Barrier degree increase");
	lemon::Arg<lemon::ArgT::Value, bry_int_t> subd = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("subd").flag('s').description("Barrier subdivision");
	lemon::Arg<lemon::ArgT::Value, bry_int_t> max_constraints = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("max-constraints").flag('c').defaultValue(10000).description("Maximum number of constraints in the adaptive problem");
	lemon::Arg<lemon::ArgT::Value, bry_int_t> n_sols = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("nsols").defaultValue(10).description("Number of solutions of the search algorithm should compare against");
	lemon::Arg<lemon::ArgT::Value, bry_int_t> iters = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("iters").defaultValue(1).description("Number of iterations to run the adaptive algorithm for");
	//lemon::Arg<lemon::ArgT::Value, bry_int_t> ada_max_subdiv = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("ada-max-subdiv").defaultValue(2l).description("Max number of sets divided each iteration (ONLY FOR ADAPTIVE)");
	lemon::Arg<lemon::ArgT::Value, bry_int_t> time_steps = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("ts").flag('t').defaultValue(10l).description("Number of time steps");
	lemon::Arg<lemon::ArgT::Value, bry_float_t> boundary_width = parser.addDef<lemon::ArgT::Value, bry_float_t>().key("boundary-width").defaultValue(0.2).description("Width of the boundary region buffer (unsafe)");
    parser.enableHelp();


    std::shared_ptr<AdaptiveProblem<DIM>> prob(new AdaptiveProblem<DIM>());

    if (dynamics_type.value() == "to_origin") {
        prob->dynamics.reset(new PolynomialDynamics<DIM>(makeUniformArray<bry_int_t, DIM>(1)));
        PolynomialDynamics<DIM>& dynamics = *prob->dynamics;
        for (bry_int_t d = 0; d < DIM; ++d) {
            std::array<bry_int_t, DIM> coeff_idx = makeUniformArray<bry_int_t, DIM>(0);
            coeff_idx[d] = 1;
            dynamics[d].coeff(coeff_idx) = 0.5;
            //std::array<bry_int_t, DIM> y_coeff_idx = makeUniformArray<bry_int_t, DIM>(0);
        }

        //TODO
        //dynamics[0].coeff(0, 0) = 0.0;
        //dynamics[0].coeff(1, 0) = 0.5;
        //dynamics[0].coeff(0, 1) = 0.0;
        //dynamics[0].coeff(1, 1) = 0.0;
        //dynamics[1].coeff(0, 0) = 0.0;
        //dynamics[1].coeff(1, 0) = 0.0;
        //dynamics[1].coeff(0, 1) = 0.5;
        //dynamics[1].coeff(1, 1) = 0.0;
    } else if (dynamics_type.value() == "random") {
        prob->dynamics.reset(new PolynomialDynamics<DIM>(makeUniformArray<bry_int_t, DIM>(dynamics_deg.value())));
        prob->dynamics->setRandom();
    }
    if (verbose) {
        INFO("Dynamics: \n" << *prob->dynamics);
    }
    
    Covariance<DIM> cov = Covariance<DIM>::Zero();
    for (bry_int_t i = 0; i < DIM; ++i) {
        cov(i, i) = 0.01;
    }
    prob->noise.reset(new AdditiveGaussianNoise<DIM>(cov));


    auto printSetBounds = [](const HyperRectangle<DIM>& set) {
        INFO("  Lower: " << set.lower_bounds.transpose());
        INFO("  Upper: " << set.upper_bounds.transpose());
    };


    HyperRectangle<DIM> workspace(-0.5, 0.5);
    workspace.lower_bounds(0) = -1.0;
    workspace.lower_bounds(1) = -0.5;
    workspace.upper_bounds(0) = 0.5;
    workspace.upper_bounds(1) = 0.5;
    prob->setWorkspace(workspace);
    if (verbose) {
        INFO("Workspace set:");
        printSetBounds(workspace);
    }

    // Init set
    // Set all higher dimensions to btw -0.1 and 0.1
    HyperRectangle<DIM> init_set(-0.01, 0.01);
    // Edit the specific values for the 2D plane we're working in
    init_set.lower_bounds(0) = -0.8;
    init_set.upper_bounds(0) = -0.6;
    init_set.lower_bounds(1) = -0.2;
    init_set.upper_bounds(1) = 0.0;
    prob->sets.insert({ConstraintType::Init, init_set});
    if (verbose) {
        INFO("Initial set:");
        printSetBounds(init_set);
    }
    
    std::list<HyperRectangle<DIM>> boundary_sets = makeRectBoundary(workspace, boundary_width.value());
    prob->insertSets(ConstraintType::Unsafe, std::make_move_iterator(boundary_sets.begin()), std::make_move_iterator(boundary_sets.end()));

    if (non_convex) {
        // Non convex unsafe regions
        HyperRectangle<DIM> upper_region(-1.0, 1.0);
        upper_region.lower_bounds(0) = -0.57;
        upper_region.upper_bounds(0) = -0.53;
        upper_region.lower_bounds(1) = -0.17;
        upper_region.upper_bounds(1) = -0.13;
        prob->sets.insert({ConstraintType::Unsafe, upper_region});
        //printSetBounds(upper_region);
        HyperRectangle<DIM> lower_region(-1.0, 1.0);
        lower_region.lower_bounds(0) = -0.57;
        lower_region.upper_bounds(0) = -0.53;
        lower_region.lower_bounds(1) = 0.28;
        lower_region.upper_bounds(1) = 0.32;
        prob->sets.insert({ConstraintType::Unsafe, lower_region});
        //printSetBounds(lower_region);
    }
    if (verbose) {
        INFO("Unsafe sets:");
        for (auto[it, end] = prob->getSets(ConstraintType::Unsafe); it != end; ++it) {
            printSetBounds(it->second);
        }
    }
    

    // Safe set

    if (!non_convex) {
        //HyperRectangle<DIM> safe_set(-1.0, 1.0);
        //safe_set.lower_bounds(0) = -1.0;
        //safe_set.upper_bounds(0) = 0.5;
        //safe_set.lower_bounds(1) = -0.5;
        //safe_set.upper_bounds(1) = 0.5;
        //prob->safe_sets.push_back(safe_set);
        prob->sets.insert({ConstraintType::Safe, workspace});
    } else {
        {
            HyperRectangle<DIM> set(-1.0, 1.0);
            set.lower_bounds(0) = -1.0;
            set.upper_bounds(0) = -0.57;
            set.lower_bounds(1) = -0.5;
            set.upper_bounds(1) = 0.5;
            prob->sets.insert({ConstraintType::Safe, set});
        }

        {
            HyperRectangle<DIM> set(-1.0, 1.0);
            set.lower_bounds(0) = -0.57;
            set.upper_bounds(0) = -0.53;
            set.lower_bounds(1) = -0.5;
            set.upper_bounds(1) = -0.17;
            prob->sets.insert({ConstraintType::Safe, set});
        }

        {
            HyperRectangle<DIM> set(-1.0, 1.0);
            set.lower_bounds(0) = -0.57;
            set.upper_bounds(0) = -0.53;
            set.lower_bounds(1) = -0.13;
            set.upper_bounds(1) = 0.28;
            prob->sets.insert({ConstraintType::Safe, set});
        }

        {
            HyperRectangle<DIM> set(-1.0, 1.0);
            set.lower_bounds(0) = -0.57;
            set.upper_bounds(0) = -0.53;
            set.lower_bounds(1) = 0.32;
            set.upper_bounds(1) = 0.5;
            prob->sets.insert({ConstraintType::Safe, set});
        }

        {
            HyperRectangle<DIM> set(-1.0, 1.0);
            set.lower_bounds(0) = -0.53;
            set.upper_bounds(0) = 0.5;
            set.lower_bounds(1) = -0.5;
            set.upper_bounds(1) = 0.5;
            prob->sets.insert({ConstraintType::Safe, set});
        }
    }
    if (verbose) {
        INFO("Safe sets:");
        for (auto[it, end] = prob->getSets(ConstraintType::Safe); it != end; ++it) {
            printSetBounds(it->second);
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

    std::shared_ptr<PolyDynamicsProblem<DIM>> prior_prob = std::make_shared<PolyDynamicsProblem<DIM>>(*prob);
    if (subd) {
        INFO("Subdividing in " << subd.value());
        prior_prob->subdivide(subd.value());
        //
        prob->subdivide(subd.value());
        //
    }

    if (export_matrices) {
        ConstraintMatrices<DIM> constraints = prob->getConstraintMatrices();
        INFO("Exporting constraint matrices...");
        writeMatrixToFile(constraints.A, "A.txt");
        writeMatrixToFile(constraints.b, "b.txt");
        INFO("Done!");
    }

    LPSolver solver(solver_id.value());

    INFO("Solving prior...");
    Timer t("total_time");
    SynthesisResult<DIM> prior_result = synthesize(solver, *prior_prob);
    INFO("Done! Solving adapted problem...");
    //prob->existing_result = &prior_result;
    prob->max_constraints = max_constraints.value();
    prob->max_ideal_solutions_found = n_sols.value();
    prob->actions = makeSubdivisonActions<DIM>();

    SynthesisResult<DIM> result = prior_result;
    prob->existing_result = &prior_result;
    //bry_float_t best_result_psafe = prior_result.p_safe;
    for (bry_int_t iter = 0; iter < iters.value(); ++iter) {
        Timer t_i("iter_time");
        INFO("Iteration " << iter + 1 << "/" << iters.value());
        result = synthesize(solver, *prob);
        INFO("Done! (time: " << t_i.now(BRY::TimeUnit::s) << ")");
        INFO("Probability of safety: " << result.p_safe);
        prob->existing_result = &result;
    }

    //if (adaptive) {
    //    result = synthesizeAdaptive(*prob, ada_iters.value(), ada_max_subdiv.value(), solver_id.value());
    //} else {
    //    result = synthesize(*prob, solver_id.value());
    //}
    double total_time = t.now(BRY::TimeUnit::s);
    INFO("Done! (total time: " << total_time << ")");
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
 
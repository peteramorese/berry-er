#include "berry/Polynomial.h"
#include "berry/MultiIndex.h"
#include "berry/Operations.h"
#include "berry/BernsteinTransform.h"
#include "lemon/ArgParser.h"
#include "HyperRectangle.h"

#include <Eigen/Dense>

#include <iomanip>
#include <random>

using namespace BRY;

int main(int argc, char** argv) {

    constexpr std::size_t DIM = 2;

	lemon::ArgParser parser(argc, argv);
	lemon::Arg<lemon::ArgT::Value, bry_int_t> barrier_deg = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("deg").flag('d').description("Barrier degree").required();
	lemon::Arg<lemon::ArgT::Value, bry_int_t> deg_increase = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("deg-inc").flag('i').defaultValue(0l).description("Barrier degree increase");
	lemon::Arg<lemon::ArgT::Value, bry_int_t> subd = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("subd").flag('s').defaultValue(0l).description("Barrier subdivision");
	lemon::Arg<lemon::ArgT::Value, bry_int_t> trials = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("trials").flag('t').defaultValue(1l).description("Number of random trials");
	lemon::Arg<lemon::ArgT::Value, bry_int_t> seed = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("seed").description("Seed the RNG");
	lemon::Arg<lemon::ArgT::Value, bry_int_t> fidelity = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("fidelity").flag('f').defaultValue(5000l).description("Set the fidelity of the empirical LB");
    parser.enableHelp();

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> real_dist(-10.0, 10.0);
    if (seed) {
        gen.seed(seed.value());
    }

    for (bry_int_t t = 0; t < trials.value(); ++t) {
        Eigen::Tensor<bry_float_t, DIM> tensor(barrier_deg.value() + 1, barrier_deg.value() + 1);
        //tensor.setRandom();
        for (std::size_t i = 0; i < tensor.size(); ++i) {
            *(tensor.data() + i) = real_dist(gen);
        }

        BRY::Polynomial<DIM> p(std::move(tensor));
        
        bry_float_t emp_lb = 1.0e300;
        bry_int_t n_emp = fidelity.value();

        //bry_float_t corner = 1.0 / (subd.get() + 1);
        bry_float_t corner = 1.0;

        for (bry_int_t i = 0; i <= n_emp; ++i) {
            for (bry_int_t j = 0; j <= n_emp; ++j) {
                bry_float_t x1 = 0.0 + i * (corner / n_emp);
                bry_float_t x2 = 0.0 + j * (corner / n_emp);
                bry_float_t val = p(x1, x2);
                if (val < emp_lb) {
                    emp_lb = val;
                }
            }
        }
        auto printSetBounds = [](const HyperRectangle<DIM>& set) {
            DEBUG("Set bounds: [" 
                << set.lower_bounds(0) << ", " 
                << set.upper_bounds(0) << ", "
                << set.lower_bounds(1) << ", " 
                << set.upper_bounds(1) << "]");
        };

        HyperRectangle<DIM> unit_set;
        unit_set.lower_bounds = {0, 0};
        unit_set.upper_bounds = corner * Eigen::Vector<bry_float_t, DIM>::Constant(1.0);

        bry_float_t min_bern_lb = 1.0e300;
        bry_float_t min_bern_ub = 1.0e300;
        std::vector<HyperRectangle<DIM>> subd_sets = unit_set.subdivide(subd.value());
        bool all_vertex = true;
        for (auto set : subd_sets) {
        //for (uint32_t i = 0; i < subd_sets.size(); ++i) {
        //    const HyperRectangle<DIM>& set = subd_sets[i];

            printSetBounds(set);
            Matrix subd_tf = set.transformationMatrix(barrier_deg.value());
            BRY::Polynomial<DIM> p_subd = transform<DIM, Basis::Power>(p, subd_tf);


            Matrix tf = BernsteinBasisTransform<DIM>::pwrToBernMatrix(p_subd.degree(), deg_increase.value());
            auto p_tf = transform<DIM, Basis::Power, Basis::Bernstein>(p_subd, tf);
            auto[bern_lb, vertex_cond] = BernsteinBasisTransform<DIM>::infBound(p_tf);
            if (!vertex_cond) {
                all_vertex = false;
            }
            bry_float_t bern_gap = BernsteinBasisTransform<DIM>::infBoundGap(p_subd, vertex_cond, deg_increase.value());
            if (bern_lb < min_bern_lb)
                min_bern_lb = bern_lb;
            if (bern_lb + bern_gap < min_bern_ub)
                min_bern_ub = bern_lb + bern_gap;
            DEBUG("Gap: [" << bern_lb << ", " << bern_lb + bern_gap << "]");
        }

        INFO("Empirical lower bound: " << emp_lb);
        INFO("Bern lower bound: [" << min_bern_lb << ", " << min_bern_ub << "]  gap size: " << min_bern_ub - min_bern_lb);
        //INFO("Bern lower bound: [" << bern_lb << ", " << bern_lb + bern_gap << "]  gap size: " << bern_gap);
        if (all_vertex) {
            INFO("Vertex condition holds");
            continue;
        }
        if (min_bern_lb > (emp_lb + 0.00001) || (min_bern_ub - 0.00001) < emp_lb) {
            ERROR("Bern gap does not hold!");
            PAUSE;
        }
    }
}

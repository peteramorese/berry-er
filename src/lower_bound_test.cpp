#include "berry/Polynomial.h"
#include "berry/MultiIndex.h"
#include "berry/Operations.h"
#include "berry/BernsteinTransform.h"
#include "ArgParser.h"
#include "HyperRectangle.h"

#include <Eigen/Dense>

#include <iomanip>
#include <random>

using namespace BRY;

int main(int argc, char** argv) {

    constexpr std::size_t DIM = 2;

	ArgParser parser(argc, argv);
	auto barrier_deg = parser.parse<bry_int_t>("deg", 'd', 3l, "Barrier degree");
	auto deg_increase = parser.parse<bry_int_t>("deg-inc", 'i', 0l, "Barrier degree increase");
	auto subd = parser.parse<bry_int_t>("subd", 's', 0l, "Barrier subdivision");
	auto trials = parser.parse<bry_int_t>("trials", 't', 1l, "Number of random trials");
	auto seed = parser.parse<bry_int_t>("seed", "Seed the RNG");
	auto fidelity = parser.parse<bry_int_t>("fidelity", 'f', 5000l, "Set the fidelity of the empirical LB");
    parser.enableHelp();

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> real_dist(-10.0, 10.0);
    if (seed.has()) {
        gen.seed(seed.get());
    }

    for (bry_int_t t = 0; t < trials.get(); ++t) {
        Eigen::Tensor<bry_float_t, DIM> tensor(barrier_deg.get() + 1, barrier_deg.get() + 1);
        tensor.setRandom();
        for (std::size_t i = 0; i < tensor.size(); ++i) {
            *(tensor.data() + i) = real_dist(gen);
        }

        BRY::Polynomial<DIM> p(std::move(tensor));
        
        bry_float_t emp_lb = 1.0e300;
        bry_int_t n_emp = fidelity.get();

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

        HyperRectangle<DIM> unit_set;
        unit_set.lower_bounds = {0, 0};
        unit_set.upper_bounds = corner * Eigen::Vector<bry_float_t, DIM>::Constant(1.0);

        bry_float_t min_bern_lb = 1.0e300;
        bry_float_t min_bern_ub = 1.0e300;
        std::vector<HyperRectangle<DIM>> subd_sets = unit_set.subdivide(subd.get());
        for (auto set : subd_sets) {

            Matrix subd_tf = set.transformationMatrix(barrier_deg.get());
            BRY::Polynomial<DIM> p_subd = transform<DIM, Basis::Power>(p, subd_tf);


            Matrix tf = BernsteinBasisTransform<DIM>::pwrToBernMatrix(p_subd.degree(), deg_increase.get());
            auto p_tf = transform<DIM, Basis::Power, Basis::Bernstein>(p_subd, tf);
            bry_float_t bern_lb = BernsteinBasisTransform<DIM>::infBound(p_tf);
            bry_float_t bern_gap = BernsteinBasisTransform<DIM>::infBoundGap(p_subd, deg_increase.get());
            if (bern_lb < min_bern_lb)
                min_bern_lb = bern_lb;
            if (bern_lb + bern_gap < min_bern_ub)
                min_bern_ub = bern_lb + bern_gap;
        }

        INFO("Empirical lower bound: " << emp_lb);
        INFO("Bern lower bound: [" << min_bern_lb << ", " << min_bern_ub << "]  gap size: " << min_bern_ub - min_bern_lb);
        //INFO("Bern lower bound: [" << bern_lb << ", " << bern_lb + bern_gap << "]  gap size: " << bern_gap);
        if (!(min_bern_lb <= (emp_lb + 0.0000001) && (min_bern_ub - 0.0000001) >= emp_lb)) {
            ERROR("Bern gap does not hold!");
            PAUSE;
        }
    }
}

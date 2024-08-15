#include "berry/Polynomial.h"
#include "berry/MultiIndex.h"
#include "berry/Operations.h"
#include "berry/BernsteinTransform.h"
#include "lemon/ArgParser.h"
#include "lemon/Random.h"
#include "HyperRectangle.h"
#include "Synthesis.h"

#include <Eigen/Dense>

#include <iomanip>
#include <random>

using namespace BRY;

int main(int argc, char** argv) {

    constexpr std::size_t DIM = 2;

	lemon::ArgParser parser(argc, argv);
	lemon::Arg<lemon::ArgT::Value, bry_int_t> barrier_deg = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("deg").flag('d').description("Barrier degree").required();
	lemon::Arg<lemon::ArgT::Value, bry_int_t> deg_increase = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("deg-inc").flag('i').defaultValue(0l).description("Barrier degree increase");
	lemon::Arg<lemon::ArgT::Value, bry_int_t> subd = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("subd").flag('s').defaultValue(2l).description("Barrier subdivision");
	lemon::Arg<lemon::ArgT::Value, bry_int_t> iters = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("iters").description("Iterations").required();
	lemon::Arg<lemon::ArgT::Value, bry_int_t> seed = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("seed").description("Seed the RNG");
	lemon::Arg<lemon::ArgT::Value, bry_int_t> fidelity = parser.addDef<lemon::ArgT::Value, bry_int_t>().key("fidelity").flag('f').defaultValue(5000l).description("Set the fidelity of the empirical LB");
	lemon::Arg<lemon::ArgT::Check> split = parser.addDef<lemon::ArgT::Check>().key("split").description("Split the sets instead of subdividing");
	lemon::Arg<lemon::ArgT::Check> export_poly = parser.addDef<lemon::ArgT::Check>().key("export-poly").flag('e').description("Export the polynomial being used (for visualization)");
	lemon::Arg<lemon::ArgT::Check> input = parser.addDef<lemon::ArgT::Check>().key("input").description("Input the split dimension");
    parser.enableHelp();

    //std::random_device rd;
    //std::mt19937_64 gen(rd());
    //std::uniform_real_distribution<double> real_dist(-10.0, 10.0);
    //if (seed) {
    //    gen.seed(seed.value());
    //}

    Eigen::Tensor<bry_float_t, DIM> tensor(barrier_deg.value() + 1, barrier_deg.value() + 1);
    Eigen::VectorXd coeff_vec(tensor.size());
    //tensor.setRandom();
    if (seed) {
        lemon::RNG::seed(seed.value());
    }
    for (std::size_t i = 0; i < tensor.size(); ++i) {
        bry_float_t random_coeff = lemon::RNG::srandd(-100.0, 100.0);
        *(tensor.data() + i) = -random_coeff;
        coeff_vec[i] = -random_coeff;
    }
    BRY::Polynomial<DIM> p(std::move(tensor));
    if (export_poly) {
        INFO("Exporting polynomial: " << p);
        writeMatrixToFile(coeff_vec, "certificate_coeffs.txt");
    }
        
    std::vector<BRY::Polynomial<DIM>> derivatives;
    derivatives.reserve(DIM);
    for (bry_int_t d = 0; d < DIM; ++d) {
        derivatives.push_back(p.derivative(d));
    }

    // Calculate the empirical lower bound
    Eigen::Vector<bry_float_t, DIM> emp_min_pt;
    bry_float_t emp_lb = 1.0e300;
    bry_int_t n_emp = fidelity.value();

    bry_float_t corner = 1.0;

    for (bry_int_t i = 0; i <= n_emp; ++i) {
        for (bry_int_t j = 0; j <= n_emp; ++j) {
            bry_float_t x1 = 0.0 + i * (corner / n_emp);
            bry_float_t x2 = 0.0 + j * (corner / n_emp);
            bry_float_t val = p(x1, x2);
            if (val < emp_lb) {
                emp_lb = val;
                emp_min_pt[0] = x1;
                emp_min_pt[1] = x2;
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

    std::list<HyperRectangle<DIM>> sets = {HyperRectangle<DIM>()};
    for (bry_int_t i = 0; i < iters.value(); ++i) {
        
        bry_float_t min_bern_lb = 1.0e300;
        bry_float_t max_gap = -1.0e300;
        //bry_float_t min_bern_ub = 1.0e300;
        //std::vector<HyperRectangle<DIM>> subd_sets = unit_set.subdivide(subd.value());
        //bool all_vertex = true;
        Eigen::Vector<bry_float_t, DIM> min_pt;
        std::array<bry_int_t, DIM> min_coeff_idx;

        auto it_to_subd = sets.end();
        bool is_min_on_vertex = false;

        NEW_LINE;
        INFO("SETS:");
        for (const auto& set : sets) {
            printSetBounds(set);
        }
        for (auto it = sets.begin(); it != sets.end(); ++it) {
            Matrix subd_tf = it->transformationMatrix(barrier_deg.value());
            BRY::Polynomial<DIM> p_subd = transform<DIM, Basis::Power>(p, subd_tf);


            Matrix phi = BernsteinBasisTransform<DIM>::pwrToBernMatrix(p_subd.degree(), deg_increase.value());
            auto p_tf = transform<DIM, Basis::Power, Basis::Bernstein>(p_subd, phi);

            std::array<bry_int_t, DIM> coeff_idx;
            auto[bern_lb, vertex_cond] = BernsteinBasisTransform<DIM>::infBound(p_tf, coeff_idx);
            //if (!vertex_cond) {
            //    all_vertex = false;
            //}
            if (vertex_cond)
                DEBUG("VERTEX COND");
            Eigen::Vector<bry_float_t, DIM> coeff_pt = BernsteinBasisTransform<DIM>::ctrlPtOnUnitBox(coeff_idx, p_tf.degree());

            bry_float_t bern_gap = BernsteinBasisTransform<DIM>::infBoundGap(p_subd, vertex_cond, deg_increase.value());
            if (bern_lb < min_bern_lb) {
                min_bern_lb = bern_lb;
                it_to_subd = it;
                is_min_on_vertex = vertex_cond;
                min_pt = coeff_pt;
                min_coeff_idx = coeff_idx;
            }
            if (bern_gap > max_gap) {
                max_gap = bern_gap;
            }
            //DEBUG("coeff_idx: (" << coeff_idx[0] << ", " << coeff_idx[1] << ")" << "   min point: [" << coeff_pt.transpose() << "]");
            INFO("LB: " << bern_lb << ", Gap: [" << bern_lb << ", " << bern_lb + bern_gap << "]");
        }
        INFO("Iter " << i << " with " << sets.size() << " sets. Lower bound: " << min_bern_lb << ", min point: " << min_pt.transpose() << "    DIFFERENCE: " << emp_lb - min_bern_lb);
        if (is_min_on_vertex) {
            WARN("Min found on vertex");
            break;
        }
        if (split) {
            INFO("Splitting set:");
            printSetBounds(*it_to_subd);
            INFO(" at [" << min_pt.transpose() << "]");
            
            // Find the dimension with smallest second derivative 
            bry_float_t min_slope = 1.0e300;
            bry_int_t split_dim;
            for (bry_int_t d = 0; d < DIM; ++d) {
                if (min_coeff_idx[d] == 0 || min_coeff_idx[d] == p.degree()) {
                    continue;
                }
                bry_float_t slope_d = std::abs(derivatives[d](min_pt));
                DEBUG("slope of dim " << d << " = " << slope_d);
                if (slope_d < min_slope) {
                    min_slope = slope_d;
                    split_dim = d;
                }
            }
            bry_int_t input_split_dim = -1;
            if (input) {
                std::cout << "Enter split dim: ";
                std::cin >> input_split_dim;
            }
            if (input_split_dim >= 0) {
                if (input_split_dim >= DIM) {
                    ERROR("Split dim out of range");
                } else {
                    split_dim = input_split_dim;
                }
            }
            DEBUG("     Splitting dim " << split_dim);
            auto[new_set_1, new_set_2] = it_to_subd->splitByPercent(split_dim, min_pt[split_dim]);
            sets.erase(it_to_subd);
            sets.push_back(new_set_1);
            sets.push_back(new_set_2);
        } else {
            INFO("Subdividing set:");
            printSetBounds(*it_to_subd);
            std::vector<HyperRectangle<DIM>> subd_sets = it_to_subd->subdivide(subd.value());
            sets.erase(it_to_subd);
            sets.insert(sets.end(), subd_sets.begin(), subd_sets.end());
        }
    }
    INFO("Empirical lower bound: " << emp_lb << " at [" << emp_min_pt.transpose() << "]");
}

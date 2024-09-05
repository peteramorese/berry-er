#include <Eigen/Dense>

#include "MatMul.h"
#include "Time.h"
#include "lemon/Logging.h"
#include "lemon/ArgParser.h"

#include <iostream>
#include <fstream>
#include <iomanip>

int main(int argc, char** argv) {

#ifdef BRYR_USE_CUBLAS

    INFO("Calibrating GPU accelerated matrix multiplication");

    lemon::ArgParser parser(argc, argv);
    lemon::Arg<lemon::ArgT::Check> verbose = parser.addDef<lemon::ArgT::Check>().flag('v').key("Verbose").description("Show all benchmarks");
    lemon::Arg<lemon::ArgT::Check> validate = parser.addDef<lemon::ArgT::Check>().key("validate").description("Check the equivalence of the matrices to ensure that cuBLAS is working");
    lemon::Arg<lemon::ArgT::Value, int> step = parser.addDef<lemon::ArgT::Value, int>().flag('s').key("step").description("Step in dimension").defaultValue(1);

    parser.enableHelp();

    INFO("Running initial dummy product on GPU...");
    {
        Eigen::MatrixXd A(10, 10);
        Eigen::MatrixXd B(10, 10);
        A.setRandom();
        B.setRandom();
        BRY::multiplyCUDA(A, B);
    }
    INFO("Done! Benchmarking...");

    int sz_thr;
    for (int dimension = 10; dimension < 10000; dimension += step.value()) {

        sz_thr = dimension * dimension;

        Eigen::MatrixXd A(dimension, dimension);
        Eigen::MatrixXd B(dimension, dimension);
        A.setRandom();
        B.setRandom();

        double time_reg, time_gpu;
        Eigen::MatrixXd C_reg, C_gpu;
        // Timer for regular matrix multiplication
        {
            BRY::Timer t("t_reg");
            C_reg = A * B;
            time_reg = t.now(BRY::TimeUnit::us);
        }

        // Timer for GPU matrix multiplication
        {   
            BRY::Timer t("t_gpu");
            C_gpu = BRY::multiplyCUDA(A, B);
            time_gpu = t.now(BRY::TimeUnit::us);
        }   

        if (validate) {
            bool eq = C_reg.isApprox(C_gpu);
            if (!eq) {
                ERROR("Matrix products are not equal, ensure cuBLAS is working");
                return 1;
            }
        }

        if (verbose) {
            INFO("Size: " << std::setw(8) << dimension * dimension << " | BLAS (us): " << std::setw(10) << time_reg << " | cuBLAS (us): " << std::setw(10) << time_gpu);
        } else {
            INFO_SMLN("Size: " << std::setw(8) << dimension * dimension << " | BLAS (us): " << std::setw(10) << time_reg << " | cuBLAS (us): " << std::setw(10) << time_gpu);
        } 
        if (time_gpu < time_reg) {
            NEW_LINE;
            INFO("Size threshold found: " << dimension * dimension);
            std::ofstream ofs("size_threshold", std::ios::trunc);

            if (ofs.is_open()) {
                ofs << sz_thr;
                ofs.close();
            } else {
                ERROR("Error writing to 'size_threshold' file");
            }
            return 0;
        }
    }
    
    WARN("Size threshold was not found between " << 100 << " and " << 10000 * 10000);
    return 1;
#else
    ERROR("BRYR_USE_CUBLAS not defined, did not build with cuBLAS capability");
    return 1;
#endif
}

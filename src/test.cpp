//#include "berry/Operations.h"
//#include "HyperRectangle.h"

//#define EIGEN_USE_BLAS
#include <Eigen/Dense>

#include "MatMul.h"
#include <iostream>

//using namespace BRY;
int main(int argc, char** argv) {
    //Eigen::setNbThreads(16);
    //Eigen::MatrixXd A = Eigen::MatrixXd::Random(10000, 50000);
    //Eigen::MatrixXd B = Eigen::MatrixXd::Random(50000, 10000);
    Eigen::MatrixXd A(4, 3);
	A << 1, 2, 3,
		4, 5, 6,
		7, 8, 9,
		10, 11, 12;
    Eigen::MatrixXd B(3, 3);
	B << 9, 8, 7,
		6, 5 ,4,
		3, 2, 1;
    //Eigen::MatrixXd C = A * B;
    Eigen::MatrixXd C = BRY_MATMUL(A, B);
	std::cout << C;
    
    //HyperRectangle<2> h1(.5, 1.2, 3);
    //HyperRectangle<2> h2(.1, 1.2, 4);

    //bry_int_t b_deg = 3;
    //DiagDegFilter<2> f(b_deg);
    //DEBUG("Matrix w/o filter: \n" << h1.transformationMatrix(b_deg));
    //NEW_LINE;
    //DEBUG("Matrix w/ filter: \n" << h1.transformationMatrix(b_deg, &f));
    //NEW_LINE;
    //Matrix new_mat = f.applyToCoeffMatrixCols(h1.transformationMatrix(b_deg));
    ////new_mat = f.applyToCoeffMatrixRows(new_mat);
    //DEBUG("matrix w/ manual filter\n" << new_mat);
    ////HyperRectangle<2> h2(.5000000000009, 1.2, 3);

    ////DEBUG("h1 < h2? " << (h1 < h2));
    ////DEBUG("h2 < h1? " << (h2 < h1));
    return 0;
}

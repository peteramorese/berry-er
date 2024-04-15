

#include <iostream>

#include <Eigen/Dense>

int main() {
    // Define a matrix
    Eigen::MatrixXd matrix(3, 3);
    matrix << 1, 2, 3,
              4, 5, 6,
              7, 8, 9;

    // Define a vector
    Eigen::VectorXd vector(3);
    vector << 10, 20, 30;

    // Set the second column of the matrix equal to the vector
    matrix.col(1) = vector;

    // Print the modified matrix
    std::cout << "Modified Matrix:\n" << matrix << std::endl;

    return 0;
}
 
#include <coin/CoinMessageHandler.hpp>
#include <coin/CoinPackedMatrix.hpp>
#include <coin/CoinModel.hpp>
#include <coin/ClpSimplex.hpp>

#include <iostream>

int main() {
    // Define the objective coefficients
    double objCoefficients[] = {3.0, 2.0};

    // Define the constraint matrix
    int constraintIndices[] = {0, 1, 0, 1};
    int constraintStarts[] = {0, 2, 4};
    double constraintCoefficients[] = {2.0, 1.0, 4.0, -5.0};
    CoinPackedMatrix matrix(false, 2, 2, 4, constraintCoefficients, constraintIndices, constraintStarts);

    // Define the variable lower bounds
    double lb[] = {0.0, 0.0};

    // Define the variable upper bounds
    double ub[] = {COIN_DBL_MAX, COIN_DBL_MAX};

    // Define the constraint lower bounds
    double constraintLowerBounds[] = {-COIN_DBL_MAX, -10.0};

    // Define the constraint upper bounds
    double constraintUpperBounds[] = {20.0, COIN_DBL_MAX};

    // Create the Coin LP solver
    ClpSimplex solver;

    // Load the problem into the solver
    solver.loadProblem(matrix, lb, ub, objCoefficients, constraintLowerBounds, constraintUpperBounds);

    // Solve the LP problem
    solver.primal();

    // Print the solution status
    std::cout << "Solution Status: ";
    switch (solver.status()) {
        case 0:
            std::cout << "Optimal" << std::endl;
            break;
        case 1:
            std::cout << "Primal Infeasible" << std::endl;
            break;
        case 2:
            std::cout << "Dual Infeasible" << std::endl;
            break;
        case 5:
            std::cout << "Stopped" << std::endl;
            break;
        default:
            std::cout << "Unknown" << std::endl;
    }

    // Print the optimal variable values
    std::cout << "Optimal Values:" << std::endl;
    const double* primalValues = solver.primalColumnSolution();
    for (int i = 0; i < 2; ++i) {
        std::cout << "x" << i << " = " << primalValues[i] << std::endl;
    }

    // Print the optimal objective function value
    std::cout << "Optimal Objective Function Value (Max): " << solver.objectiveValue() << std::endl;

    return 0;
}
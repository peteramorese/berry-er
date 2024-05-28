#include <iostream>
#include <memory>

#include "ortools/linear_solver/linear_solver.h"

namespace operations_research {
void LinearProgrammingExample() {
  std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("SCIP"));
  if (!solver) {
    LOG(WARNING) << "SCIP solver unavailable.";
    return;
  }

  const double infinity = solver->infinity();
  // x and y are non-negative variables.
  MPVariable* const x = solver->MakeNumVar(0.0, infinity, "x");
  MPVariable* const y = solver->MakeNumVar(0.0, infinity, "y");
  LOG(INFO) << "Number of variables = " << solver->NumVariables();

  // x + 2*y <= 14.
  MPConstraint* const c0 = solver->MakeRowConstraint(-infinity, 14.0);
  c0->SetCoefficient(x, 1);
  c0->SetCoefficient(y, 2);

  // 3*x - y >= 0.
  MPConstraint* const c1 = solver->MakeRowConstraint(0.0, infinity);
  c1->SetCoefficient(x, 3);
  c1->SetCoefficient(y, -1);

  // x - y <= 2.
  MPConstraint* const c2 = solver->MakeRowConstraint(-infinity, 2.0);
  c2->SetCoefficient(x, 1);
  c2->SetCoefficient(y, -1);
  LOG(INFO) << "Number of constraints = " << solver->NumConstraints();

  // Objective function: 3x + 4y.
  MPObjective* const objective = solver->MutableObjective();
  objective->SetCoefficient(x, 3);
  objective->SetCoefficient(y, 4);
  objective->SetMaximization();

  const MPSolver::ResultStatus result_status = solver->Solve();
  // Check that the problem has an optimal solution.
  if (result_status != MPSolver::OPTIMAL) {
    LOG(FATAL) << "The problem does not have an optimal solution!";
  }

  LOG(INFO) << "Solution:";
  LOG(INFO) << "Optimal objective value = " << objective->Value();
  LOG(INFO) << x->name() << " = " << x->solution_value();
  LOG(INFO) << y->name() << " = " << y->solution_value();
}
}  // namespace operations_research

int main(int argc, char** argv) {
  operations_research::LinearProgrammingExample();
  return EXIT_SUCCESS;
}
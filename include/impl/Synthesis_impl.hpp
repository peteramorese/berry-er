#pragma once

#include "Synthesis.h"

#include "berry/BernsteinTransform.h"

#include <fstream>
#include <iomanip>

template <std::size_t DIM>
void BRY::SynthesisResult<DIM>::removeFilter() {
    if (!filter) {
        return;
    }

    bry_int_t n_coeffs = pow(barrier_deg + 1, DIM);
    ASSERT(b_values.size() != n_coeffs, "Supplied coefficient vector is already in square-degree form (no filter is applied)");

    Vector sq_coeffs = Vector::Zero(n_coeffs);
    bry_int_t i = 0;

    for (auto col_midx = mIdxW(DIM, barrier_deg + 1); !col_midx.last(); ++col_midx) {
        if (!filter->remove(col_midx.begin())) {
            sq_coeffs(col_midx.inc().wrappedIdx()) = b_values(i++);
        }
    }

    b_values = sq_coeffs;
}

template <std::size_t DIM>
BRY::SynthesisResult<DIM> BRY::synthesize(const PolyDynamicsProblem<DIM>& problem, const std::string& solver_id) {
    auto constraints = problem.getConstraintMatrices();
    return synthesize(constraints, problem.time_horizon, solver_id);
}

template <std::size_t DIM>
BRY::SynthesisResult<DIM> BRY::synthesize(const ConstraintMatrices<DIM>& constraints, bry_int_t time_horizon, const std::string& solver_id) {
    LPSolver solver(solver_id);
    solver.setConstraintMatrices(constraints.A, constraints.b);

    BRY::SynthesisResult<DIM> result(constraints.barrier_deg, constraints.filter);
    static_cast<LPSolver::Result&>(result) = solver.solve(time_horizon);
    return result;
}

template <std::size_t DIM>
BRY::SynthesisResult<DIM> BRY::synthesizeAdaptive(PolyDynamicsProblem<DIM> problem, bry_int_t max_iter, bry_int_t subdiv_per_iter, const std::string& solver_id) {
    ASSERT(subdiv_per_iter >= 1, "subdiv_per_iter must be geq than 1");

    bry_int_t iter = 0;
    while (true) {
        LPSolver solver(solver_id);
        ConstraintMatrices<DIM> constraints = problem.getConstraintMatrices(true);
        solver.setConstraintMatrices(constraints.A, constraints.b);

        BRY::SynthesisResult<DIM> result(constraints.barrier_deg, constraints.filter);
        static_cast<LPSolver::Result&>(result) = solver.solve(problem.time_horizon);

        NEW_LINE;
        INFO("[Iteration " << iter + 1 << " / " << max_iter << "] p_safe = " << result.p_safe << " (time: " << result.comp_time <<")");

        if (++iter == max_iter) {
            return result;
        }

        if (result.isFilterApplied()) {
            result.removeFilter();
        }
        Polynomial<DIM> barrier(result.b_values);


        // Compute the smallest max_subdiv_per_iter robustness values
        Vector soln_vec = solver.getSolnVector();
        Vector robustness_values = constraints.computeRobustnessVec(soln_vec);
        ASSERT(subdiv_per_iter < problem.numSets(), "Number of subdivisions per iter is more than the number of sets");
        std::vector<bry_int_t> robustness_indices(robustness_values.size());
        for (bry_int_t i = 0; i < robustness_indices.size(); ++i) {
            robustness_indices[i] = i;
        }

        // Sort the robustness values
        std::sort(robustness_indices.begin(), robustness_indices.end(),
            [&robustness_values](bry_int_t lhs, bry_int_t rhs){return robustness_values[lhs] < robustness_values[rhs];});

        // Determine the sets that have the lowest robustness values
        auto comp = [](const typename std::list<HyperRectangle<DIM>>::iterator& lhs, const typename std::list<HyperRectangle<DIM>>::iterator& rhs) -> bool {
            return &(*lhs) < &(*rhs); 
        };
        std::set<typename std::list<HyperRectangle<DIM>>::iterator, decltype(comp)> unique_iterators(comp);

        struct UniqueSetToSubd {
            typename std::list<HyperRectangle<DIM>>::iterator it; // It pointing to the set to subd
            ConstraintType type; // List to edit the element in
            bry_float_t bound_gap; // Conservativeness of the bound in this set
        };
        std::vector<UniqueSetToSubd> bound_gap_sorted_iters;
        bound_gap_sorted_iters.reserve(subdiv_per_iter);

        for (bry_int_t robustness_idx : robustness_indices) {
            const ConstraintID& id = constraints.constraint_ids[robustness_idx];
            auto it = problem.lookupSetFromConstraint(id);
            bool inserted = unique_iterators.insert(it).second;
            if (inserted) {
                Matrix tf = constraints.transformation_matrices->at(id);
                Polynomial<DIM> constrained_barrier = transform(barrier, tf);
                bry_float_t bound_gap = BernsteinBasisTransform<DIM>::infBoundGap(constrained_barrier, false, problem.degree_increase);
                bound_gap_sorted_iters.push_back(UniqueSetToSubd{it, id.type, bound_gap});
                if ((bound_gap_sorted_iters.size() >= subdiv_per_iter) && (std::abs(robustness_values[robustness_idx]) > BRY_FLOAT_DIFF_TOL)) {
                    break;
                }
            }
        }

        // Sort the iters by decreasing bound gap to subd the ones with largest gap first
        std::sort(bound_gap_sorted_iters.begin(), bound_gap_sorted_iters.end(), 
            [](const UniqueSetToSubd& lhs, const UniqueSetToSubd& rhs){return lhs.bound_gap > rhs.bound_gap;});

        // Subdivide the marked sets by removing the original set and adding the subdivided sets
        auto subdivideAndReplace = [](std::list<HyperRectangle<DIM>>& list_of_sets, std::list<HyperRectangle<DIM>>::iterator it_to_subd) {
            // Subdivide in 2
            std::vector<HyperRectangle<DIM>> subdivisions = it_to_subd->subdivide(2);
            // Erase the original set
            auto following_it = list_of_sets.erase(it_to_subd);
            // Insert the subdivided sets
            list_of_sets.insert(following_it, subdivisions.begin(), subdivisions.end());
        };

        for (const UniqueSetToSubd& set_to_subd : bound_gap_sorted_iters) {
            switch (set_to_subd.type) {
                case ConstraintType::Workspace:
                    subdivideAndReplace(problem.workspace_sets, set_to_subd.it);
                    break;
                case ConstraintType::Init:
                    subdivideAndReplace(problem.init_sets, set_to_subd.it);
                    break;
                case ConstraintType::Unsafe:
                    subdivideAndReplace(problem.unsafe_sets, set_to_subd.it);
                    break;
                case ConstraintType::Safe:
                    subdivideAndReplace(problem.safe_sets, set_to_subd.it);
            }
        }
    }
}

void BRY::writeMatrixToFile(const Matrix& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        // Write matrix dimensions
        file << matrix.rows() << " " << matrix.cols() << std::endl;

        // Write matrix data
        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = 0; j < matrix.cols(); ++j) {
                file << std::fixed << std::setprecision(40) << matrix(i, j);
                if (j < matrix.cols() - 1) file << " ";
            }
            file << std::endl;
        }

        file.close();
    } else {
        std::cerr << "Unable to open file";
    }
}
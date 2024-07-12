#pragma once

#include "Synthesis.h"

#include "berry/BernsteinTransform.h"

#include <fstream>
#include <iomanip>

template <std::size_t DIM>
BRY::ConstraintMatrices<DIM>::ConstraintMatrices(bry_int_t n_constraints, bry_int_t n_vars, bry_int_t barrier_deg_)
    : A(n_constraints, n_vars)
    , b(n_constraints)
    , barrier_deg(barrier_deg_)
    , constraint_ids(n_constraints)
    , m_diag_deg(false)
{
    ASSERT(A.cols() == pow(barrier_deg + 1, DIM) + 2, "Number of vars does not match barrier degree + 2");
}

template <std::size_t DIM>
void BRY::ConstraintMatrices<DIM>::toDiagonalDegree() {
    if (m_diag_deg) {
        WARN("Already in using diagonal degree form");
        return;
    }
    bry_int_t n_coeffs = pow(barrier_deg + 1, DIM);

    auto removeCol = [](Eigen::MatrixXd& m, int col) {
        Eigen::MatrixXd new_m(m.rows(), m.cols() - 1);
        new_m << m.leftCols(col), m.rightCols(m.cols() - col - 1);
        m = new_m;
    };

    bry_int_t cols_removed = 0;

    for (auto col_midx = mIdxW(DIM, barrier_deg + 1); !col_midx.last(); ++col_midx) {
        if (std::accumulate(col_midx.begin(), col_midx.end(), 0) > barrier_deg) {
            removeCol(A, col_midx.inc().wrappedIdx() - cols_removed++);
        }
    }
    m_diag_deg = true;
}

template <std::size_t DIM>
BRY::Vector BRY::ConstraintMatrices<DIM>::computeRobustnessVec(const Vector& soln_vec) const {
    return A * soln_vec - b;
}

template <std::size_t DIM>
void BRY::PolyDynamicsProblem<DIM>::setWorkspace(const HyperRectangle<DIM>& workspace) {
    workspace_sets = {workspace};
}

template <std::size_t DIM>
void BRY::PolyDynamicsProblem<DIM>::subdivide(uint32_t subdivision) {
    if (subdivision < 2) {
        WARN("Subdivision is less than 2 (no effect)");
        return;
    }
    auto divide = [&] (std::list<HyperRectangle<DIM>>& subd_sets) {
        const std::list<HyperRectangle<DIM>> original_sets = subd_sets;
        subd_sets.clear();
        for (const auto& set : original_sets) {
            std::vector<HyperRectangle<DIM>> subd_sets_to_insert = set.subdivide(subdivision);
            subd_sets.insert(subd_sets.end(), subd_sets_to_insert.begin(), subd_sets_to_insert.end());
        }
    };
    divide(workspace_sets);
    divide(init_sets);
    divide(safe_sets);
    divide(unsafe_sets);
}

template <std::size_t DIM>
BRY::bry_int_t BRY::PolyDynamicsProblem<DIM>::numSets() const {
    return workspace_sets.size() + init_sets.size() + unsafe_sets.size() + safe_sets.size();
}

template <std::size_t DIM>
const BRY::ConstraintMatrices<DIM> BRY::PolyDynamicsProblem<DIM>::getConstraintMatrices() const {

    Matrix Phi_m = BernsteinBasisTransform<DIM>::pwrToBernMatrix(barrier_deg, degree_increase);

    bry_int_t n_cols = Phi_m.cols() + 2;

    // Workspace
    if (workspace_sets.empty())
        WARN("No workspace was provided");
    Matrix ws_coeffs(workspace_sets.size() * Phi_m.rows(), n_cols);
    Vector ws_lower_bound = Vector::Zero(ws_coeffs.rows());
    bry_int_t i = 0;
    for (const HyperRectangle<DIM>& set : workspace_sets) {
        Matrix beta_coeffs = Phi_m * set.transformationMatrix(barrier_deg);
        Matrix coeffs(beta_coeffs.rows(), n_cols);
        coeffs << beta_coeffs, Vector::Zero(beta_coeffs.rows()), Vector::Zero(beta_coeffs.rows());
        ws_coeffs.block(Phi_m.rows() * i++, 0, Phi_m.rows(), n_cols) = coeffs;
    }

    // Initial sets
    if (init_sets.empty())
        WARN("No initial sets were provided");
    Matrix init_coeffs(init_sets.size() * Phi_m.rows(), n_cols);
    Vector init_lower_bound = Vector::Zero(init_coeffs.rows());
    i = 0;
    for (const HyperRectangle<DIM>& set : init_sets) {
        Matrix beta_coeffs = -Phi_m * set.transformationMatrix(barrier_deg);
        Matrix coeffs(beta_coeffs.rows(), n_cols);
        coeffs << beta_coeffs, Vector::Ones(beta_coeffs.rows()), Vector::Zero(beta_coeffs.rows());
        init_coeffs.block(Phi_m.rows() * i++, 0, Phi_m.rows(), n_cols) = coeffs;
    }

    // Unsafe sets
    if (unsafe_sets.empty())
        WARN("No unsafe sets were provided");
    Matrix unsafe_coeffs(unsafe_sets.size() * Phi_m.rows(), n_cols);
    Vector unsafe_lower_bound = Vector::Ones(unsafe_coeffs.rows());
    i = 0;
    for (const HyperRectangle<DIM>& set : unsafe_sets) {
        Matrix beta_coeffs = Phi_m * set.transformationMatrix(barrier_deg);
        Matrix coeffs(beta_coeffs.rows(), n_cols);
        coeffs << beta_coeffs, Vector::Zero(beta_coeffs.rows()), Vector::Zero(beta_coeffs.rows());
        unsafe_coeffs.block(Phi_m.rows() * i++, 0, Phi_m.rows(), n_cols) = coeffs;
    }

    // Safe sets
    if (safe_sets.empty())
        WARN("No safe sets were provided");
    //DEBUG("Dynamics power matrix: \n" << dynamics->dynamicsPowerMatrix(barrier_deg));
    //DEBUG("Noise matrix: \n" << noise->additiveNoiseMatrix(barrier_deg));
    Matrix F_expec_Gamma = dynamics->dynamicsPowerMatrix(barrier_deg) * noise->additiveNoiseMatrix(barrier_deg);
    bry_int_t p = dynamics->composedDegree(barrier_deg);
    Matrix Phi_p = BernsteinBasisTransform<DIM>::pwrToBernMatrix(p, degree_increase);
    Vector gamma_coeffs = Vector::Ones(Phi_p.cols());
    ASSERT(F_expec_Gamma.rows() == Phi_p.cols(), "Dimension mismatch between F and Phi (p)");

    Matrix safe_coeffs(safe_sets.size() * Phi_p.rows(), n_cols);
    Vector safe_lower_bound = Vector::Zero(safe_coeffs.rows());
    i = 0;

    Matrix deg_lift_tf = makeDegreeChangeTransform<DIM>(barrier_deg, p);

    for (const HyperRectangle<DIM>& set : safe_sets) {
        //DEBUG("Set transformation:\n" << set.transformationMatrix(p));
        Matrix beta_coeffs = 
            -Phi_p * set.transformationMatrix(p) * (F_expec_Gamma - deg_lift_tf);
        
        Matrix coeffs(beta_coeffs.rows(), n_cols);
        coeffs << beta_coeffs, Vector::Zero(beta_coeffs.rows()), Vector::Ones(beta_coeffs.rows());
        safe_coeffs.block(Phi_p.rows() * i++, 0, Phi_p.rows(), n_cols) = coeffs;
    }

    BRY::ConstraintMatrices<DIM> constraint_matrices(ws_coeffs.rows() + init_coeffs.rows() + unsafe_coeffs.rows() + safe_coeffs.rows(), n_cols, barrier_deg);

    constraint_matrices.A << ws_coeffs, init_coeffs, unsafe_coeffs, safe_coeffs;
    constraint_matrices.b << ws_lower_bound, init_lower_bound, unsafe_lower_bound, safe_lower_bound;

    // Fill the constraint IDs
    auto it = constraint_matrices.constraint_ids.begin();
    for (bry_int_t set_i = 0; set_i < workspace_sets.size(); ++set_i) {
        std::fill(it, std::next(it, Phi_m.rows()), ConstraintID{ConstraintType::Workspace, set_i});
        std::advance(it, Phi_m.rows());
    }
    for (bry_int_t set_i = 0; set_i < init_sets.size(); ++set_i) {
        std::fill(it, std::next(it, Phi_m.rows()), ConstraintID{ConstraintType::Init, set_i});
        std::advance(it, Phi_m.rows());
    }
    for (bry_int_t set_i = 0; set_i < unsafe_sets.size(); ++set_i) {
        std::fill(it, std::next(it, Phi_m.rows()), ConstraintID{ConstraintType::Unsafe, set_i});
        std::advance(it, Phi_m.rows());
    }
    for (bry_int_t set_i = 0; set_i < safe_sets.size(); ++set_i) {
        std::fill(it, std::next(it, Phi_p.rows()), ConstraintID{ConstraintType::Safe, set_i}); // Advance by Phi_p rows instead of Phi_m
        std::advance(it, Phi_p.rows());
    }
    
    if (diag_deg) {
        constraint_matrices.toDiagonalDegree();
    }
    return constraint_matrices;
}

template <std::size_t DIM>
std::list<BRY::HyperRectangle<DIM>>::iterator BRY::PolyDynamicsProblem<DIM>::lookupSetFromConstraint(const ConstraintID& id) {
    switch (id.type) {
        case ConstraintType::Workspace: {
            #ifdef BRY_ENABLE_BOUNDS_CHECK
                ASSERT(id.set_idx < workspace_sets.size(), "Set idx out of bounds (workspace sets)");
            #endif
            return std::next(workspace_sets.begin(), id.set_idx);
        }
        case ConstraintType::Init: {
            #ifdef BRY_ENABLE_BOUNDS_CHECK
                ASSERT(id.set_idx < init_sets.size(), "Set idx out of bounds (init sets)");
            #endif
            return std::next(init_sets.begin(), id.set_idx);
        }
        case ConstraintType::Unsafe: {
            #ifdef BRY_ENABLE_BOUNDS_CHECK
                ASSERT(id.set_idx < unsafe_sets.size(), "Set idx out of bounds (unsafe sets)");
            #endif
            return std::next(unsafe_sets.begin(), id.set_idx);
        }
        case ConstraintType::Safe: {
            #ifdef BRY_ENABLE_BOUNDS_CHECK
                ASSERT(id.set_idx < safe_sets.size(), "Set idx out of bounds (safe sets)");
            #endif
            return std::next(safe_sets.begin(), id.set_idx);
        }
    }
    throw std::invalid_argument("ID is invalid");
}

template <std::size_t DIM>
void BRY::SynthesisResult<DIM>::fromDiagonalDegree() {
    bry_int_t n_coeffs = pow(barrier_deg + 1, DIM);
    ASSERT(b_values.size() != n_coeffs, "Supplied coefficient vector is already in square-degree form");
    Vector sq_coeffs = Vector::Zero(n_coeffs);

    bry_int_t i = 0;

    for (auto col_midx = mIdxW(DIM, barrier_deg + 1); !col_midx.last(); ++col_midx) {
        if (std::accumulate(col_midx.begin(), col_midx.end(), 0) <= barrier_deg) {
            sq_coeffs(col_midx.inc().wrappedIdx()) = b_values(i++);
        }
    }

    b_values = sq_coeffs;
}

template <std::size_t DIM>
BRY::SynthesisResult<DIM> BRY::synthesize(const PolyDynamicsProblem<DIM>& problem, const std::string& solver_id) {
    LPSolver solver(solver_id);
    auto constraints = problem.getConstraintMatrices();
    return synthesize(constraints, problem.time_horizon, solver_id);
}

template <std::size_t DIM>
BRY::SynthesisResult<DIM> BRY::synthesize(const ConstraintMatrices<DIM>& constraints, bry_int_t time_horizon, const std::string& solver_id) {
    LPSolver solver(solver_id);
    solver.setConstraintMatrices(constraints.A, constraints.b);

    BRY::SynthesisResult<DIM> result(constraints.diagDeg(), constraints.barrier_deg);
    static_cast<LPSolver::Result&>(result) = solver.solve(time_horizon);
    return result;
}

template <std::size_t DIM>
BRY::SynthesisResult<DIM> BRY::synthesizeAdaptive(PolyDynamicsProblem<DIM> problem, bry_int_t max_iter, bry_int_t subdiv_per_iter, const std::string& solver_id) {
    ASSERT(subdiv_per_iter >= 1, "subdiv_per_iter must be geq than 1");
    BRY::SynthesisResult<DIM> result(problem.diag_deg, problem.barrier_deg);
    for (bry_int_t iter = 0; iter < max_iter; ++iter) {
        LPSolver solver(solver_id);
        ConstraintMatrices<DIM> constraints = problem.getConstraintMatrices();
        solver.setConstraintMatrices(constraints.A, constraints.b);

        static_cast<LPSolver::Result&>(result) = solver.solve(problem.time_horizon);

        NEW_LINE;
        INFO("Iteration " << iter + 1 << " / " << max_iter << " p_safe = " << result.p_safe);

        if (iter == max_iter - 1) {
            break;
        }

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
        std::set<typename std::list<HyperRectangle<DIM>>::iterator, decltype(comp)> ws_sets_to_subd(comp);
        std::set<typename std::list<HyperRectangle<DIM>>::iterator, decltype(comp)> init_sets_to_subd(comp);
        std::set<typename std::list<HyperRectangle<DIM>>::iterator, decltype(comp)> unsf_sets_to_subd(comp);
        std::set<typename std::list<HyperRectangle<DIM>>::iterator, decltype(comp)> sf_sets_to_subd(comp);
        bry_int_t num_sets_to_subd = 0;
        for (bry_int_t robustness_idx : robustness_indices) {
            const ConstraintID& id = constraints.constraint_ids[robustness_idx];
            auto it = problem.lookupSetFromConstraint(id);
            bool inserted = false;
            switch (id.type) {
                case ConstraintType::Workspace:
                    inserted = ws_sets_to_subd.insert(it).second;
                    DEBUG("Subd ws set: " <<  id.set_idx);
                    break;
                case ConstraintType::Init: 
                    inserted = init_sets_to_subd.insert(it).second;
                    DEBUG("Subd init set: " <<  id.set_idx);
                    break;
                case ConstraintType::Unsafe:
                    inserted = unsf_sets_to_subd.insert(it).second;
                    DEBUG("Subd unsafe set: " <<  id.set_idx);
                    break;
                case ConstraintType::Safe:
                    inserted = sf_sets_to_subd.insert(it).second;
                    DEBUG("Subd safe set: " <<  id.set_idx);
                    break;
            }

            if (inserted) {
                ++num_sets_to_subd;
                if (num_sets_to_subd >= subdiv_per_iter) {
                    break;
                }
            }
        }

        // Subdivide the marked sets by removing the original set and adding the subdivided sets
        auto subdivideAndReplace = [](std::list<HyperRectangle<DIM>>& list_of_sets, std::list<HyperRectangle<DIM>>::iterator it_to_subd) {
            //DEBUG("removing it: " << &*it_to_subd);
            //bool found = false;
            //for (auto it = list_of_sets.begin(); it != list_of_sets.end(); ++it) {
            //    if (it == it_to_subd) {
            //        found = true;
            //        break;
            //    }
            //}
            //ASSERT(found, "Iterator not found");

            // Subdivide in 2
            std::vector<HyperRectangle<DIM>> subdivisions = it_to_subd->subdivide(2);
            // Erase the original set
            auto following_it = list_of_sets.erase(it_to_subd);
            // Insert the subdivided sets
            list_of_sets.insert(following_it, subdivisions.begin(), subdivisions.end());
        };

        if (iter >= 1) {
            //for (uint32_t i = 0; i < robustness_indices.size(); ++i) {
            //    ConstraintID id = constraints.constraint_ids[robustness_indices[i]];
            //    DEBUG("Robustness val" << i << " " << robustness_values[robustness_indices[i]] << " set type: " << (uint32_t)id.type << " idx: " << id.set_idx);
            //}

            DEBUG("TESTING! ACTUALLY DIVIDING type: " << (uint32_t)ConstraintType::Safe << " set: 2");
            auto it = problem.lookupSetFromConstraint(ConstraintID{ConstraintType::Safe, 2});
            subdivideAndReplace(problem.safe_sets, it);

            //DEBUG("TESTING! ACTUALLY DIVIDING Unsafe set 0");
            //auto it = problem.lookupSetFromConstraint(ConstraintID{ConstraintType::Unsafe, 0});
            //subdivideAndReplace(problem.unsafe_sets, it);

            continue;
        }

        for (auto it : ws_sets_to_subd) {
            subdivideAndReplace(problem.workspace_sets, it);
        }
        for (auto it : init_sets_to_subd) {
            subdivideAndReplace(problem.init_sets, it);
        }
        for (auto it : unsf_sets_to_subd) {
            subdivideAndReplace(problem.unsafe_sets, it);
        }
        for (auto it : sf_sets_to_subd) {
            subdivideAndReplace(problem.safe_sets, it);
        }
        //auto printSetBounds = [](const HyperRectangle<DIM>& set) {
        //    DEBUG("Set bounds: [" 
        //        << set.lower_bounds(0) << ", " 
        //        << set.upper_bounds(0) << ", "
        //        << set.lower_bounds(1) << ", " 
        //        << set.upper_bounds(1) << "]");
        //};

        //DEBUG("Original list:");
        //for (auto set : problem.unsafe_sets) {
        //    printSetBounds(set);
        //}
        //NEW_LINE;
        //auto test_it = ++problem.unsafe_sets.begin();
        //DEBUG("Set to remove:");
        //printSetBounds(*test_it);
        //subdivideAndReplace(problem.unsafe_sets, test_it);
        //NEW_LINE;

        //DEBUG("Edited list:");
        //for (auto set : problem.unsafe_sets) {
        //    printSetBounds(set);
        //}

    }
    return result;
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
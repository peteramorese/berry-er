#pragma once

#include "Problem.h"

#include <set>

namespace BRY {

template <std::size_t DIM>
class AdaptiveProblem : public PolyDynamicsProblem<DIM> {
    public:
        /// @brief Array of virtual actions (owning) pointers that the algorithm will use. Actions will be deleted on destruction
        std::vector<Action<DIM>*> actions;

        /// @brief Maximum number of constraints in the problem
        bry_int_t max_constraints;

        /// @brief Maximum number of solutions the search algorithms encounters before returning the best one.
        bry_int_t max_ideal_solutions_found = 10;
        
        /// @brief If not null, adapt the subdivision to the existing barrier up to max_constraints
        const LPSolver::Result* existing_result = nullptr;

    public:
        virtual const ConstraintMatrices<DIM> getConstraintMatrices() const override;

};

}

#include "impl/AdaptiveProblem_impl.hpp"
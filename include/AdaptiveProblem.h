#pragma once

#include "Problem.h"

namespace BRY {

template <std::size_t DIM>
class AdaptiveProblem : public PolyDynamicsProblem<DIM> {
    public:
        /// @brief Maximum number of constraints in the problem
        bry_int_t max_constraints;
        
        /// @brief If not null, adapt the subdivision to the existing barrier up to max_constraints
        const LPSolver::Result* existing_result = nullptr;

    public:
        virtual const ConstraintMatrices<DIM> getConstraintMatrices(bool store_tf_matrices = false) const override;
};

}

#include "impl/AdaptiveProblem_impl.hpp"
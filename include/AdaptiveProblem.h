#pragma once

#include "Problem.h"

namespace BRY {

template <std::size_t DIM>
struct Action {
    virtual apply(SetDefinitions<DIM>& set_defs, typename SetDefinitions<DIM>::Iterator application_set) const = 0;
};

/* Actions that can be used when constructing the adaptive problem */

/// @brief Divide along a given dimension
/// @tparam DIM Dimension of the problem
/// @tparam DIV_DIM Dimension in which the set should be divided
template <std::size_t DIM, std::size_t DIV_DIM>
struct Divide : public Action {
    virtual apply(SetDefinitions<DIM>& set_defs, typename SetDefinitions<DIM>::Iterator application_set) const override;
};

/// @brief Increase the bernstein degree of a set
/// @tparam DIM Dimension of the problem
template <std::size_t DIM>
struct IncreaseDegree : public Action {
    virtual apply(SetDefinitions<DIM>& set_defs, typename SetDefinitions<DIM>::Iterator application_set) const override;

    /// @brief Amount to increase the degree by each time the action is applied
    bry_int_t increase = 1;
};

template <std::size_t DIM>
class AdaptiveProblem : public PolyDynamicsProblem<DIM> {
    public:
        /// @brief Array of virtual actions (owning) pointers that the algorithm will use. Actions will be deleted on destruction
        std::vector<Action*> actions;

        /// @brief Maximum number of constraints in the problem
        bry_int_t max_constraints;
        
        /// @brief If not null, adapt the subdivision to the existing barrier up to max_constraints
        const LPSolver::Result* existing_result = nullptr;

    public:
        virtual const ConstraintMatrices<DIM> getConstraintMatrices() const override;

    private:
        struct AppliedAction {
            /// @brief Enumeration of the set the action was applied to
            bry_int_t set_enumeration;
            /// @brief Pointer to the action pointer inside `actions`
            Action** action;
        };

        struct State {
            /// @brief Ordered sequence of action layers such that all actions 
            std::list<std::vector<AppliedAction>> action_layers;
        };

        //struct State : public SetDefinitions<DIM> {
        //    //std::map<typename SetDefinitions<DIM>::ConstIterator, bry_float_t> min_robustness_value;

        //    public:

        //};
};

}

#include "impl/AdaptiveProblem_impl.hpp"
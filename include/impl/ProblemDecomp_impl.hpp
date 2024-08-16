#pragma once

#include "ProblemDecompSearch.h"

template <std::size_t DIM>
ProblemDecompSearch<DIM>::ProblemDecompSearch(const AdaptiveProblem* problem)
    : m_problem(problem)
{}

template <std::size_t DIM>
void BRY::ProblemDecompSearch<DIM>::reset() {
    m_unique_sets.clear();
    m_unique_states.clear();
    m_expansion_set.clear();
    m_solution_states_encountered = 0;
    m_solution_state = nullptr;
}

template <std::size_t DIM>
BRY::ProblemDecompSearch<DIM>::State::State(const State& other) 
    : sets(other.sets)
    , min_robustness(other.min_robustness)
    , n_total_constraints(other.n_total_constraints)
{
    min_robustness_set = sets.find(*other.min_robustness_set);
}

template <std::size_t DIM>
const BRY::ProblemDecompSearch<DIM>::State* BRY::ProblemDecompSearch<DIM>::search() {
    // Convert the given problem into an initial state
    State init_state;
    init_state.min_robustness = 1e100;
    init_state.n_total_constraints = 0;
    for (Set set : this->sets) {
        // Insert the set into the set registry
        auto[qset_ptr, unq_inserted] = insertUniqueSet(std::move(set));
        auto[it, init_state_inserted] = init_state.insert(qset_ptr);

        ASSERT(!(unq_inserted && !init_state_inserted), "Set was inserted into registry, but initial search state already had it");

        // If the set was inserted into the unique container, it has never been encountered before, so calculate robustness
        if (unq_inserted) {
            calculateRobustness(*qset_ptr);
        } else {
            WARN("Duplicate set found in problem defintion (Set type: " << set.first << ")");
        }

        if (qset_ptr->min_robustness < init_state.min_robustness) {
            init_state.min_robustness = qset_ptr->min_robustness;
        }
        init_state.n_total_constraints += q_set_ptr->n_constraints;
    }

    // Terminate of the initial state already exceeds the max constraints
    if (init_state.n_total_constraints > m_probem->max_constraints) {
        WARN("Original set definitions exceed max constraints (search terminating)");
        return nullptr;
    }
}

template <std::size_t DIM>
std::pair<typename BRY::ProblemDecompSearch<DIM>::QuantifiedSet*, bool> BRY::ProblemDecompSearch<DIM>::insertUniqueSet(Set&& set) {
    // Insert with zero robustness and constraints for now, if the insertion takes place, the correct values will be calculated after
    auto[it, inserted] = m_unique_sets.insert({std::move(set), 0.0, 0});
    return std::make_pair(&*it, inserted);
}

template <std::size_t DIM>
void BRY::ProblemDecompSearch<DIM>::expandState(const State* curr_state, const std::vector<Action<DIM>*>& actions) {
    for (Action<DIM>* action : actions) {
        // Copy the state
        State new_state = *curr_state;

        // Apply the action to the set and get the new sets
        std::vector<HyperRectangle<DIM>> new_sets;
        action->apply((*new_state.min_robustness_set)->set, new_sets);

        // Subtract the constraints that the set was contributing 
        new_state.n_total_constraints -= *new_state.min_robustness_set->n_constraints;

        // Erase the old set from the new state, since we are replacing it with new_sets
        new_state.erase(new_state.min_robustness_set);

        // Add each new set to the unique_sets and new state
        for (HyperRectangle<DIM>& new_set : new_sets) {
            // Try inserting the set into the registry of unique sets
            auto[new_qset_ptr, unq_inserted] = insertUniqueSet(std::move(new_set));
            // Insert the new set into the new state if its not a duplicate
            auto[it, new_state_inserted] = new_state.insert(new_qset_ptr);

            ASSERT(!(unq_inserted && !new_state_inserted), "Set was inserted into registry, but a state already had it");

            // If the set was inserted into the unique container, it has never been encountered before, so calculate robustness
            if (unq_inserted) {
                calculateRobustness(*new_qset_ptr);
            }

            // If the set was inserted into the state, we need to add the number of constraints the added set is contributing
            if (new_state_inserted) {
                new_state.n_total_constraints += new_qset_ptr->n_constraints;
            }
        }
        
        // If the constraint cap has been breached, propose the current state as a solution candidate.
        if (new_state.n_total_constraints >= max_constraints) {
            proposeSolutionState(curr_state);
            continue;
        }

        // Try inserting the state, if it has already been seen, this will return false
        auto[unq_state_it, inserted] = m_unique_states.insert(std::move(new_state));

        // If the state has not been seen, then we need to calculate the new robustness values and add it to expansion set
        if (inserted) {
            // Find the min robustness element
            auto comp = [] (const QuantifiedSet* lhs, const QuantifiedSet* rhs) {return lhs->min_robustness < rhs->min_robustness;};
            unq_state_it->min_robustness_set = std::min_element(sets.begin(), sets.end(), comp);
            // Reset the min robustness value to the robustness of the found element
            unq_state_it->min_robustness =  *unq_state_it->min_robustness_set->min_robustness;
        }
    }
}


template <std::size_t DIM>
void BRY::ProblemDecompSearch<DIM>::proposeSolutionState(const State* state) {
    if (!!m_solution_state) { // If a solution state exists
        ++m_solution_states_encountered;
        if (state->min_robustness > m_solution_state->min_robustness) {
            m_solution_state = state;
        }
    } else { // otherwise set the first solution state
        m_solution_state = state;
    }
}
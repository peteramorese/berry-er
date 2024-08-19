#pragma once

#include "AdaptiveProblem.h"

#include "berry/BernsteinTransform.h"

/* Actions */

template <std::size_t DIM, std::size_t DIV_DIM>
void BRY::Divide<DIM, DIV_DIM>::apply(const HyperRectangle<DIM>& old_set, const Eigen::Vector<bry_float_t, DIM>& normalized_split_point, std::vector<HyperRectangle<DIM>>& new_sets) const {
    std::pair<HyperRectangle<DIM>, HyperRectangle<DIM>> split_sets = old_set.splitByPercent(DIV_DIM, normalized_split_point[DIV_DIM]);
    new_sets.reserve(2);
    new_sets.push_back(std::move(split_sets.first));
    new_sets.push_back(std::move(split_sets.second));
}

template <std::size_t DIM>
void BRY::IncreaseDegree<DIM>::apply(const HyperRectangle<DIM>& old_set, const Eigen::Vector<bry_float_t, DIM>& normalized_split_point, std::vector<HyperRectangle<DIM>>& new_sets) const {
    HyperRectangle<DIM> deg_incr_set = old_set;
    deg_incr_set.bernstein_deg_incr += increase;
    new_sets = {std::move(deg_incr_set)};
}

/* AdaptiveProblem */

template <std::size_t DIM>
const BRY::ConstraintMatrices<DIM> BRY::AdaptiveProblem<DIM>::getConstraintMatrices() {
    // If there is no existing result, it must be the first iteration and thus return the normal constraint matrices
    if (!existing_result) {
        return PolyDynamicsProblem<DIM>::getConstraintMatrices();
    }
    INFO("Creating constraint matrices");

    // Degree of the composed polynomial
    m_p = this->dynamics->composedDegree(this->barrier_deg);
    // Product F * E[Gamma]
    m_F_expec_Gamma = this->dynamics->dynamicsPowerMatrix(this->barrier_deg) * this->noise->additiveNoiseMatrix(this->barrier_deg);
    // Degree lift transform for the subtraction of F * E[Gamma] - B
    m_deg_lift_tf = makeDegreeChangeTransform<DIM>(this->barrier_deg, m_p);
    // Number of optimization variables
    m_n_cols = BRY::pow(this->barrier_deg + 1, DIM) + 2;

    if (!!existing_result) {
        m_soln_vec = Vector(m_n_cols);
        m_soln_vec.segment(0, existing_result->b_values.size()) = existing_result->b_values;
        m_soln_vec[m_n_cols - 2] = existing_result->eta;
        m_soln_vec[m_n_cols - 1] = existing_result->gamma;
    }

    const State* search(); 

    //return BRY::ConstraintMatrices<DIM>();
}

template <std::size_t DIM>
void BRY::AdaptiveProblem<DIM>::reset() {
    m_unique_sets.clear();
    m_unique_states.clear();
    m_expansion_set.clear();
    m_solution_states_encountered = 0;
    m_solution_state = nullptr;
}

template <std::size_t DIM>
BRY::AdaptiveProblem<DIM>::State::State(const State& other) 
    : sets(other.sets)
    , min_robustness(other.min_robustness)
    , n_total_constraints(other.n_total_constraints)
{
    min_robustness_set = sets.find(*other.min_robustness_set);
}

template <std::size_t DIM>
const BRY::AdaptiveProblem<DIM>::State* BRY::AdaptiveProblem<DIM>::search() {
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
        init_state.n_total_constraints += qset_ptr->n_constraints;
    }

    // Terminate of the initial state already exceeds the max constraints
    if (init_state.n_total_constraints > max_constraints) {
        WARN("Original set definitions exceed max constraints (search terminating)");
        return nullptr;
    }

    while (m_solution_states_encountered < max_ideal_solutions_found) {
        // Pop the expansion state off the top
        auto top_it = m_expansion_set.begin();
        const State* expansion_state = *top_it;

        // Remove state from the expansion set
        m_expansion_set.erase(top_it);

        if (*(expansion_state->min_robustness_set)->vertex_condition) {
            WARN("Tried to on vertex condition, terminating search");
            return m_solution_state;
        }
        expandState(expansion_state, actions);
    }
    return m_solution_state;
}

template <std::size_t DIM>
std::pair<typename BRY::AdaptiveProblem<DIM>::QuantifiedSet*, bool> BRY::AdaptiveProblem<DIM>::insertUniqueSet(Set&& set) {
    // Insert with zero robustness and constraints for now, if the insertion takes place, the correct values will be calculated after
    auto[it, inserted] = m_unique_sets.insert({std::move(set), 0.0, 0});
    return std::make_pair(&*it, inserted);
}

template <std::size_t DIM>
void BRY::AdaptiveProblem<DIM>::expandState(const State* curr_state, const std::vector<Action<DIM>*>& actions) {
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
            unq_state_it->min_robustness_set = std::min_element(new_state.sets.begin(), new_state.sets.end(), comp);
            // Reset the min robustness value to the robustness of the found element
            unq_state_it->min_robustness =  *unq_state_it->min_robustness_set->min_robustness;
        }
    }
}


template <std::size_t DIM>
void BRY::AdaptiveProblem<DIM>::proposeSolutionState(const State* state) {
    if (!!m_solution_state) { // If a solution state exists
        ++m_solution_states_encountered;
        if (state->min_robustness > m_solution_state->min_robustness) {
            m_solution_state = state;
        }
    } else { // otherwise set the first solution state
        m_solution_state = state;
    }
}

template <std::size_t DIM>
void BRY::AdaptiveProblem<DIM>::calculateRobustness(QuantifiedSet& qset, Matrix* A_ptr, Vector* b_ptr) {
    Matrix A_local;
    Vector b_local;
    Matrix& A = (!!A_ptr) ? *A_ptr : A_local;
    Matrix& b = (!!b_ptr) ? *b_ptr : b_local;
    bry_float_t lower_bound = 0.0;

    ConstraintType set_type = qset.set.first;
    if (set_type != ConstraintType::Safe) {
        // Eta coeffs are 1 if the set is an initial set, otherwise they are zero
        bry_float_t eta_coeff = static_cast<bry_float_t>(set_type == ConstraintType::Init);
        // Lower bound is zero unless the set is of type unsafe, in which case the lb is 1
        lower_bound = static_cast<bry_float_t>(set_type == ConstraintType::Unsafe);
        // If the set is an initial set, then negate the b coefficients
        bry_float_t coeff_multiplier = (set_type == ConstraintType::Init) ? -1.0 : 1.0;

        Matrix tf = coeff_multiplier * qset.set.second.transformationMatrix(this->barrier_deg);
        Matrix coeffs = getPhim(qset.set.second.bernstein_deg_incr) * tf;

        A.resize(coeffs.rows(), m_n_cols);

        //            b         eta                                           gamma
        A << coeffs, Vector::Constant(coeffs.rows(), eta_coeff), Vector::Zero(coeffs.rows());
        b = Vector::Constant(A.rows(), lower_bound);
    } else {
        Matrix tf = -qset.set.second.transformationMatrix(m_p) * (m_F_expec_Gamma - m_deg_lift_tf);
        Matrix coeffs = getPhip(qset.set.second.bernstein_deg_incr) * tf;

        A.resize(coeffs.rows(), m_n_cols);

        //            b         eta                            gamma
        A << coeffs, Vector::Zero(coeffs.rows()), Vector::Ones(coeffs.rows());
        b = Vector::Zero(A.rows());
    }

    // Create a polynomial for determining the lower bound and control point index
    Polynomial<DIM, Basis::Bernstein> p(A * m_soln_vec);

    std::array<bry_int_t, DIM>& coefficient_idx;
    auto[inf_of_p, vertex_cond] = BernsteinBasisTransform<DIM>::infBound(p, coefficient_idx);

    qset.min_robustness = inf_of_p - lower_bound;
    qset.normalized_split_point = BernsteinBasisTransform<DIM>::ctrlPtOnUnitBox(coefficient_idx, p.degree());
    qset.vertex_condition = vertex_cond;
    qset.n_constraints = A.rows();
}

template <std::size_t DIM>
const BRY::Matrix& BRY::AdaptiveProblem<DIM>::getPhim(bry_int_t bernstein_deg_incr) {
    auto it = m_Phi_m.find(bernstein_deg_incr);
    if (it == m_Phi_m.end()) {
        it = m_Phi_m.emplace(bernstein_deg_incr, BernsteinBasisTransform<DIM>::pwrToBernMatrix(this->barrier_deg, bernstein_deg_incr)).first;
    } 
    return it->second;
}

template <std::size_t DIM>
const BRY::Matrix& BRY::AdaptiveProblem<DIM>::getPhip(bry_int_t bernstein_deg_incr) {
    auto it = m_Phi_p.find(bernstein_deg_incr);
    if (it == m_Phi_p.end()) {
        it = m_Phi_p.emplace(bernstein_deg_incr, BernsteinBasisTransform<DIM>::pwrToBernMatrix(m_p, bernstein_deg_incr)).first;
    } 
    return it->second;
}
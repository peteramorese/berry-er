#pragma once

#include "AdaptiveProblem.h"

#include "berry/BernsteinTransform.h"

/* Actions */

template <std::size_t DIM>
void BRY::Divide<DIM>::apply(const HyperRectangle<DIM>& old_set, const Eigen::Vector<bry_float_t, DIM>& normalized_split_point, std::vector<HyperRectangle<DIM>>& new_sets) const {
    std::pair<HyperRectangle<DIM>, HyperRectangle<DIM>> split_sets = old_set.splitByPercent(div_dim, normalized_split_point[div_dim]);
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

template <std::size_t DIM>
std::vector<BRY::Action<DIM>*> BRY::makeSubdivisonActions() {
    std::vector<BRY::Action<DIM>*> actions(DIM);
    for (bry_int_t d = 0; d < DIM; ++d) {
        actions[d] = new Divide<DIM>(d);
    }
    return actions;
}

/* AdaptiveProblem */

template <std::size_t DIM>
BRY::AdaptiveProblem<DIM>::~AdaptiveProblem() {
    for (Action<DIM>* action : actions) {
        delete action;
    }
}

template <std::size_t DIM>
const BRY::ConstraintMatrices<DIM> BRY::AdaptiveProblem<DIM>::getConstraintMatrices() {
    // If there is no existing result, it must be the first iteration and thus return the normal constraint matrices
    if (!existing_result) {
        return PolyDynamicsProblem<DIM>::getConstraintMatrices();
    }
    INFO("Creating constraint matrices using existing result");

    if (actions.empty()) {
        WARN("No set construction actions provided");
        return PolyDynamicsProblem<DIM>::getConstraintMatrices();
    }

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

    const State* ideal_state = search(); 
    INFO("Optimal state found (" << ideal_state->n_total_constraints << " constraints)");

    BRY::ConstraintMatrices<DIM> constraint_matrices(ideal_state->n_total_constraints, m_n_cols, this->barrier_deg);

    bry_int_t constraint_idx = 0;
    for (QSetIt qset_it : ideal_state->sets) {
        auto[A, b] = calculateConstraintMatrices(qset_it->first);
        constraint_matrices.A.block(constraint_idx, 0, A.rows(), A.cols()) = A;
        constraint_matrices.b.segment(constraint_idx, b.size()) = b;
        constraint_idx += A.rows();
    }

    return constraint_matrices;
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
        auto[qset_it, unq_inserted] = insertUniqueSet(std::move(set));
        auto[it, init_state_inserted] = init_state.sets.insert(qset_it);

        ASSERT(!(unq_inserted && !init_state_inserted), "Set was inserted into registry, but initial search state already had it");

        // If the set was inserted into the unique container, it has never been encountered before, so calculate robustness
        if (unq_inserted) {
            calculateRobustness(qset_it->first, qset_it->second);
        } else {
            WARN("Duplicate set found in problem defintion (Set type: " << set.first << ")");
        }

        if (qset_it->second.min_robustness < init_state.min_robustness) {
            init_state.min_robustness = qset_it->second.min_robustness;
        }
        init_state.n_total_constraints += qset_it->second.n_constraints;
    }

    // Terminate of the initial state already exceeds the max constraints
    if (init_state.n_total_constraints > max_constraints) {
        WARN("Original set definitions exceed max constraints (search terminating)");
        return nullptr;
    }

    ASSERT(m_unique_states.empty(), "State container has not been cleared before search");

    auto it = m_unique_states.insert(std::move(init_state)).first;
    m_expansion_set.insert(&*it);

    while (m_solution_states_encountered < max_ideal_solutions_found) {
        // Pop the expansion state off the top
        auto top_it = m_expansion_set.begin();
        const State* expansion_state = *top_it;

        INFO_SMLN("Robustness: " << expansion_state->min_robustness);
        // Remove state from the expansion set
        m_expansion_set.erase(top_it);

        if ((*expansion_state->min_robustness_set)->second.vertex_condition) {
            WARN("Tried to on vertex condition, terminating search");
            return m_solution_state;
        }
        expandState(expansion_state, actions);
    }
    return m_solution_state;
}

template <std::size_t DIM>
std::pair<typename BRY::AdaptiveProblem<DIM>::QSetIt, bool> BRY::AdaptiveProblem<DIM>::insertUniqueSet(Set&& set) {
    // Insert with zero robustness and constraints for now, if the insertion takes place, the correct values will be calculated after
    return m_unique_sets.insert(std::make_pair(std::move(set), SetProperties{}));
}

template <std::size_t DIM>
void BRY::AdaptiveProblem<DIM>::expandState(const State* curr_state, const std::vector<Action<DIM>*>& actions) {
    for (Action<DIM>* action : actions) {
        // Copy the state
        State new_state = *curr_state;

        // Apply the action to the set and get the new sets
        std::vector<HyperRectangle<DIM>> new_sets;
        action->apply((*new_state.min_robustness_set)->first.second, (*new_state.min_robustness_set)->second.normalized_split_point, new_sets);

        // Subtract the constraints that the set was contributing 
        new_state.n_total_constraints -= (*new_state.min_robustness_set)->second.n_constraints;

        // Pull out the constraint type from the set that is being replaced
        ConstraintType constraint_type = (*new_state.min_robustness_set)->first.first;

        // Erase the old set from the new state, since we are replacing it with new_sets
        new_state.sets.erase(new_state.min_robustness_set);

        // Add each new set to the unique_sets and new state
        for (HyperRectangle<DIM>& new_set : new_sets) {
            // Try inserting the set into the registry of unique sets
            auto[qset_it, unq_sets_inserted] = insertUniqueSet(std::make_pair(constraint_type, std::move(new_set)));

            if (unq_sets_inserted) {
                calculateRobustness(qset_it->first, qset_it->second);
            }

            // Insert the new set into the new state if its not a duplicate
            auto[it, new_state_sets_inserted] = new_state.sets.insert(qset_it);

            // If the set was inserted into the state, we need to add the number of constraints the added set is contributing
            if (new_state_sets_inserted) {
                new_state.n_total_constraints += qset_it->second.n_constraints;
            }
        }
        
        // If the constraint cap has been breached, propose the current state as a solution candidate.
        if (new_state.n_total_constraints >= max_constraints) {
            proposeSolutionState(curr_state);
            continue;
        }

        // Check if the state is new to check if we need to calculate min robustness
        bool state_is_new = !m_unique_states.contains(new_state);
        //auto[unq_state_it, inserted] = m_unique_states.insert(std::move(new_state));

        // If the state has not been seen, then we need to calculate the new robustness values and add it to expansion set
        if (state_is_new) {
            // Find the min robustness element
            auto comp = [] (const QSetIt& lhs, const QSetIt& rhs) {return lhs->second.min_robustness < rhs->second.min_robustness;};
            new_state.min_robustness_set = std::min_element(new_state.sets.begin(), new_state.sets.end(), comp);
            // Reset the min robustness value to the robustness of the found element
            new_state.min_robustness =  (*new_state.min_robustness_set)->second.min_robustness;
            auto[it, inserted] = m_unique_states.insert(std::move(new_state));
            ASSERT(inserted, "New state has not been encountered before but was not inserted");
            m_expansion_set.insert(&*it);
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
std::pair<BRY::Matrix, BRY::Vector> BRY::AdaptiveProblem<DIM>::calculateConstraintMatrices(const Set& set) {
    Matrix A;
    Vector b;
    bry_float_t lower_bound = 0.0;

    ConstraintType set_type = set.first;
    if (set_type != ConstraintType::Safe) {
        // Eta coeffs are 1 if the set is an initial set, otherwise they are zero
        bry_float_t eta_coeff = static_cast<bry_float_t>(set_type == ConstraintType::Init);
        // Lower bound is zero unless the set is of type unsafe, in which case the lb is 1
        lower_bound = static_cast<bry_float_t>(set_type == ConstraintType::Unsafe);
        // If the set is an initial set, then negate the b coefficients
        bry_float_t coeff_multiplier = (set_type == ConstraintType::Init) ? -1.0 : 1.0;

        Matrix tf = coeff_multiplier * set.second.transformationMatrix(this->barrier_deg);
        Matrix coeffs = getPhim(set.second.bernstein_deg_incr) * tf;

        A.resize(coeffs.rows(), m_n_cols);

        //            b         eta                                           gamma
        A << coeffs, Vector::Constant(coeffs.rows(), eta_coeff), Vector::Zero(coeffs.rows());
        b = Vector::Constant(A.rows(), lower_bound);
    } else {
        Matrix tf = -set.second.transformationMatrix(m_p) * (m_F_expec_Gamma - m_deg_lift_tf);
        Matrix coeffs = getPhip(set.second.bernstein_deg_incr) * tf;

        A.resize(coeffs.rows(), m_n_cols);

        //            b         eta                            gamma
        A << coeffs, Vector::Zero(coeffs.rows()), Vector::Ones(coeffs.rows());
        b = Vector::Zero(A.rows());
    }
    return std::make_pair(std::move(A), std::move(b));
}

template <std::size_t DIM>
void BRY::AdaptiveProblem<DIM>::calculateRobustness(const Set& set, SetProperties& properties) {
    auto[A, b] = calculateConstraintMatrices(set);
    bry_float_t lower_bound = b[0];

    // Create a polynomial for determining the lower bound and control point index
    Polynomial<DIM, Basis::Bernstein> p(A * m_soln_vec);

    std::array<bry_int_t, DIM> coefficient_idx;
    auto[inf_of_p, vertex_cond] = BernsteinBasisTransform<DIM>::infBound(p, coefficient_idx);

    properties.min_robustness = inf_of_p - lower_bound;
    properties.normalized_split_point = BernsteinBasisTransform<DIM>::ctrlPtOnUnitBox(coefficient_idx, p.degree());
    properties.vertex_condition = vertex_cond;
    properties.n_constraints = A.rows();
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
#pragma once

#include "Problem.h"
#include "MatMul.h"

#include "berry/BernsteinTransform.h"

#include <fstream>
#include <iomanip>


bool BRY::ConstraintID::operator<(const ConstraintID& other) const {
    if (type == other.type) {
        return set_idx < other.set_idx;
    } else {
        return type < other.type;
    }
}

template <std::size_t DIM>
bool BRY::SetDefinitions<DIM>::operator<(const SetDefinitions& other) const {
    return sets < other.sets;
}

template <std::size_t DIM>
std::pair<typename BRY::SetDefinitions<DIM>::ConstIterator, typename BRY::SetDefinitions<DIM>::ConstIterator> BRY::SetDefinitions<DIM>::getSets(ConstraintType set_type) const {
    return std::make_pair(sets.lower_bound(set_type), sets.upper_bound(set_type));
}

template <std::size_t DIM>
template <typename IT>
void BRY::SetDefinitions<DIM>::insertSets(ConstraintType set_type, IT begin, IT end) {
    for (IT it = begin; it != end; ++it) {
        sets.insert(std::make_pair(set_type, *it));
    }
}

template <std::size_t DIM>
void BRY::SetDefinitions<DIM>::setWorkspace(const HyperRectangle<DIM>& workspace) {
    sets.insert(std::make_pair(ConstraintType::Workspace, workspace));
}

template <std::size_t DIM>
void BRY::SetDefinitions<DIM>::subdivide(uint32_t subdivision) {
    if (subdivision < 2) {
        WARN("Subdivision is less than 2 (no effect)");
        return;
    }
    const std::multimap<ConstraintType, HyperRectangle<DIM>> original_sets = sets;
    sets.clear();
    for (const auto&[type, original_set] : original_sets) {
        std::vector<HyperRectangle<DIM>> subd_sets_to_insert = original_set.subdivide(subdivision);
        insertSets(type, subd_sets_to_insert.begin(), subd_sets_to_insert.end());
    }
}


template <std::size_t DIM>
BRY::ConstraintMatrices<DIM>::ConstraintMatrices(bry_int_t n_constraints, bry_int_t n_vars, bry_int_t barrier_deg_)
    : A(n_constraints, n_vars)
    , b(n_constraints)
    , barrier_deg(barrier_deg_)
    , constraint_sets(n_constraints)
{
    ASSERT(A.cols() == pow(barrier_deg + 1, DIM) + 2, "Number of vars does not match barrier degree + 2");
}

template <std::size_t DIM>
BRY::ConstraintMatrices<DIM>::ConstraintMatrices(bry_int_t n_constraints, bry_int_t n_vars, bry_int_t barrier_deg_, const std::shared_ptr<MonomialFilter<DIM>>& filter_)
    : A(n_constraints, n_vars)
    , b(n_constraints)
    , barrier_deg(barrier_deg_)
    , filter(filter_)
    , constraint_sets(n_constraints)
{
    if (!!filter) {
        ASSERT(n_vars == filter->nRemainingMonoms() + 2, "Number of vars does not filtered monomials + 2");
    } else {
        ASSERT(A.cols() == pow(barrier_deg + 1, DIM) + 2, "Number of vars does not match barrier degree + 2");
    }
}

template <std::size_t DIM>
BRY::Vector BRY::ConstraintMatrices<DIM>::computeRobustnessVec(const Vector& soln_vec) const {
    return A * soln_vec - b;
}

template <std::size_t DIM>
const BRY::ConstraintMatrices<DIM> BRY::PolyDynamicsProblem<DIM>::getConstraintMatrices() {
    INFO("Creating constraint matrices");

    initMatrixDefinitions();

    // Keep track of the number of constraints for each type
    std::array<bry_int_t, 4> n_constraints = makeUniformArray<bry_int_t, 4>(0);
    //bry_int_t n_ws_constraints = 0, n_init_constraints = 0, n_unsafe_constraints = 0, n_safe_constraints = 0;

    for (const auto&[set_type, set] : this->sets) {
        // Safe set constraints use the Phi_p container, all other constraints use Phi_m. Accumulate the number of rows
        n_constraints[set_type] += (set_type != ConstraintType::Safe) ? getPhim(set.bernstein_deg_incr).rows() : getPhip(set.bernstein_deg_incr).rows();
    }

    BRY::ConstraintMatrices<DIM> constraint_matrices(std::accumulate(n_constraints.begin(), n_constraints.end(), 0), m_n_cols, barrier_deg, filter);

    bry_int_t constraint_idx = 0;
    for (typename SetDefinitions<DIM>::ConstIterator it = this->sets.begin(); it != this->sets.end(); ++it) {
        auto[A, b] = calculateSetConstraints(it->first, it->second);
        constraint_matrices.A.block(constraint_idx, 0, A.rows(), A.cols()) = A;
        constraint_matrices.b.segment(constraint_idx, b.size()) = b;
        constraint_idx += A.rows();
    }
    return constraint_matrices;
}

template <std::size_t DIM>
void BRY::PolyDynamicsProblem<DIM>::increaseDegree(bry_int_t deg_incr) {
    for (auto&[ct, set] : this->sets) {
        set.bernstein_deg_incr = deg_incr;
    }
}

template <std::size_t DIM>
void BRY::PolyDynamicsProblem<DIM>::refineResult(LPSolver::Result& result, bry_int_t eta_iterations, bry_int_t gamma_iterations) const {
    //if (result.isFilterApplied()) {
    //    result.removeFilter();
    //}

    auto createSolnVec = [&result] (bry_float_t eta, bry_float_t gamma) -> Vector {
        Vector soln_vec(result.b_values.size() + 2);
        soln_vec.segment(0, result.b_values.size()) = result.b_values;
        soln_vec[soln_vec.size() - 2] = eta;
        soln_vec[soln_vec.size() - 1] = gamma;
        return soln_vec;
    };

    auto calculatePolynomial = [&result, this] (ConstraintType ct, const HyperRectangle<DIM>& set) -> Polynomial<DIM, Basis::Bernstein> {
        if (ct == ConstraintType::Init) {
            Vector p_vector = -this->getPhim(set.bernstein_deg_incr) * set.transformationMatrix(this->barrier_deg, filter.get()) * result.b_values;
            return Polynomial<DIM, Basis::Bernstein>(std::move(p_vector));
        } else if (ct == ConstraintType::Safe) {
            Matrix tf = set.transformationMatrix(this->m_p);
            //DEBUG("b4 prod");
            Matrix coeffs = -this->getPhip(set.bernstein_deg_incr) * tf * (m_F_expec_Gamma_minus_I);
            Vector p_vector = coeffs * result.b_values;
            //Vector p_vector = -this->getPhip(set.bernstein_deg_incr) * set.transformationMatrix(this->m_p) * (m_F_expec_Gamma_minus_I) * result.b_values;

            //DEBUG("min coeff: " << p_vector.minCoeff());
            return Polynomial<DIM, Basis::Bernstein>(std::move(p_vector));
        }
        throw std::invalid_argument("Unrecognized constraint type");
    };

    auto calculateRobustness = [&result] (ConstraintType ct, const Polynomial<DIM, Basis::Bernstein>& p) -> std::pair<bry_float_t, bool> {
        // Returns lower bound (-eta, -gamma)
        std::pair<bry_float_t, bool> inf_bound = BernsteinBasisTransform<DIM>::infBound(p);
        if (ct == ConstraintType::Init) {
            // Robustness = (0 if vertex condition) (-(-eta) - original_eta)
            inf_bound.first = result.eta + inf_bound.first;
            return inf_bound;
        } else if (ct == ConstraintType::Safe) {
            // Robustness = (0 if vertex condition) (original_gamma -(-gamma))
            //DEBUG("upper bound: " << -inf_bound.first << " gamma: " << result.gamma);
            inf_bound.first = result.gamma + inf_bound.first;
            return inf_bound;
        }
        throw std::invalid_argument("Unrecognized constraint type");
    };

    // Maps robustness to the set flagged by if it obtains vertex condition
    std::multimap<bry_float_t, std::pair<HyperRectangle<DIM>, bool>> sorted_init_sets;
    std::multimap<bry_float_t, std::pair<HyperRectangle<DIM>, bool>> sorted_safe_sets;

    // Create the initial solution vector from the quantities in the original result
    Vector soln_vec = createSolnVec(result.eta, result.gamma);

    //
    //Polynomial<DIM> p(result.b_values);
    //HyperRectangle<DIM> init_set = this->sets.find(ConstraintType::Init)->second;
    //p = transform(p, init_set.transformationMatrix(barrier_deg, filter.get()));
    //auto p_bern = transform<DIM, Basis::Power, Basis::Bernstein>(-p, BernsteinBasisTransform<DIM>::pwrToBernMatrix(barrier_deg, 100));
    //std::array<bry_int_t, DIM> coefficient_idx;
    //auto[ub, vc] = BernsteinBasisTransform<DIM>::infBound(p_bern);
    //DEBUG("max of p over init set: " << -ub << " vertex condition? " << vc);
    //PAUSE;
    //

    auto refine = [&] (std::multimap<bry_float_t, std::pair<HyperRectangle<DIM>, bool>>& sorted_container, bry_float_t max_iterations, ConstraintType constraint_type, const char* set_name) {
        // Insert the init sets into the sorted container
        //DEBUG("b4 insert");
        for (auto[it, end_it] = this->getSets(constraint_type); it != end_it; ++it) {
            //auto[A, b] = calculateSetConstraints(constraint_type, it->second);
            //bry_float_t robustness = calculateRobustness(A, b, soln_vec);

            auto[robustness, vertex_condition] = calculateRobustness(constraint_type, calculatePolynomial(constraint_type, it->second));
            sorted_container.insert(std::make_pair(robustness, std::make_pair(it->second, vertex_condition)));
            //DEBUG("robustness on insertion: " << robustness << " vertex condtion? " << vertex_condition);
        }
        //DEBUG("af insert");

        for (bry_int_t iters = 0; iters < max_iterations; ++iters) {
            // Retrieve the lowest robustness set
            auto top_it = sorted_container.begin(); 
            //DEBUG(" top_it robustness: " << top_it->first);

            // If the lowest robustness set has the vertex condition, return, can't do anything better
            if (top_it->second.second) {
                break;
            }
            
            // Subdivide it
            std::vector<HyperRectangle<DIM>> new_sets = top_it->second.first.subdivide();

            // Erase it from the container since we are inserting the subdivisions
            sorted_container.erase(top_it);
            
            for (HyperRectangle<DIM>& new_set : new_sets) {
                // Calclulate the robustness
                //auto[A, b] = calculateSetConstraints(ConstraintType::Init, new_set);
                //DEBUG("A: \n" << A << "\n b: \n" << b);
                //DEBUG("af calc new set constr");
                auto[robustness, vertex_condition] = calculateRobustness(constraint_type, calculatePolynomial(constraint_type, new_set));
                //auto p = calculatePolynomial(constraint_type, new_set);
                //for (bry_int_t i = 0; i < p.tensor().size(); ++i) {
                //    std::cout << *(p.tensor().data() + i) << " ";
                //}
                //NEW_LINE;
                
                sorted_container.insert(std::make_pair(robustness, std::make_pair(std::move(new_set), vertex_condition)));
                //bry_float_t robustness = calculateRobustness(A, b, soln_vec);
                //DEBUG("   constr new set: l" << new_set.lower_bounds[0] << ", " << new_set.lower_bounds[1] << " u" << new_set.upper_bounds[0] << ", " << new_set.upper_bounds[1]);
                //DEBUG("   robustness of new set: " << robustness << " vc? " << vertex_condition);
                //PAUSE;

                // Insert the set into the container with its robustness 
                //sorted_container.insert(std::make_pair(robustness, std::move(new_set)));
            }
            INFO_SMLN("Refining " << set_name << " sets (iteration " << std::setw(6) << iters + 1 << " / " << max_iterations << ") | Robustness: " << std::setw(15) << sorted_container.begin()->first << " | Number of sets: " << std::setw(6) << sorted_container.size());
        }
    };

    refine(sorted_init_sets, eta_iterations, ConstraintType::Init, "Init");
    refine(sorted_safe_sets, gamma_iterations, ConstraintType::Safe, "Safe");

    NEW_LINE;
    INFO("Done!");

    INFO("Eta robustness: " << sorted_init_sets.begin()->first << ", gamma robustness: " << sorted_safe_sets.begin()->first);
    bry_float_t new_eta = result.eta - sorted_init_sets.begin()->first;
    bry_float_t new_gamma = result.gamma - sorted_safe_sets.begin()->first;

    INFO("New eta: " << new_eta << ", new gamma: " << new_gamma);
    if (new_eta < result.eta) {
        result.eta = new_eta;
    }
    if (new_gamma < result.gamma) {
        result.gamma = new_gamma;
    }
}

template <std::size_t DIM>
void BRY::PolyDynamicsProblem<DIM>::initMatrixDefinitions() {
    // Degree of the composed polynomial
    m_p = dynamics->composedDegree(barrier_deg);
    // Product F * E[Gamma]
    Matrix F_expec_Gamma = dynamics->dynamicsPowerMatrix(barrier_deg) * noise->additiveNoiseMatrix(barrier_deg);
    // Subtract degree lift transform to get (F * E[Gamma] - I) * b
    m_F_expec_Gamma_minus_I = F_expec_Gamma - makeDegreeChangeTransform<DIM>(barrier_deg, m_p);

    if (!filter) {
        // Number of optimization variables
        m_n_cols = BRY::pow(barrier_deg + 1, DIM) + 2;
    } else {
        ASSERT(barrier_deg == filter->barrierDeg(), "Barrier degree used to consruct filter does not match that of the problem");
        // If there is a filter, remove the coresponding columns
        m_F_expec_Gamma_minus_I = filter->applyToCoeffMatrixCols(m_F_expec_Gamma_minus_I);
        m_n_cols = filter->nRemainingMonoms() + 2;
    }
}

template <std::size_t DIM>
const BRY::Matrix& BRY::PolyDynamicsProblem<DIM>::getPhim(bry_int_t bernstein_deg_incr) const {
    auto it = m_Phi_m.find(bernstein_deg_incr);
    if (it == m_Phi_m.end()) {
        it = m_Phi_m.emplace(bernstein_deg_incr, BernsteinBasisTransform<DIM>::pwrToBernMatrix(this->barrier_deg, bernstein_deg_incr)).first;
    } 
    return it->second;
}

template <std::size_t DIM>
const BRY::Matrix& BRY::PolyDynamicsProblem<DIM>::getPhip(bry_int_t bernstein_deg_incr) const {
    auto it = m_Phi_p.find(bernstein_deg_incr);
    if (it == m_Phi_p.end()) {
        it = m_Phi_p.emplace(bernstein_deg_incr, BernsteinBasisTransform<DIM>::pwrToBernMatrix(m_p, bernstein_deg_incr)).first;
    } 
    return it->second;
}

template <std::size_t DIM>
std::pair<BRY::Matrix, BRY::Vector> BRY::PolyDynamicsProblem<DIM>::calculateSetConstraints(ConstraintType constraint_type, const HyperRectangle<DIM>& set) const {
    Matrix A;
    Vector b;
    bry_float_t lower_bound = 0.0;
    //DEBUG("calculating set constraints for type: " << constraint_type);
    if (constraint_type != ConstraintType::Safe) {
        // Eta coeffs are 1 if the set is an initial set, otherwise they are zero
        bry_float_t eta_coeff = static_cast<bry_float_t>(constraint_type == ConstraintType::Init);
        // Lower bound is zero unless the set is of type unsafe, in which case the lb is 1
        lower_bound = static_cast<bry_float_t>(constraint_type == ConstraintType::Unsafe);
        // If the set is an initial set, then negate the b coefficients
        bry_float_t coeff_multiplier = (constraint_type == ConstraintType::Init) ? -1.0 : 1.0;

        Matrix coeffs = coeff_multiplier * BRY_MATMUL(this->getPhim(set.bernstein_deg_incr), set.transformationMatrix(this->barrier_deg, filter.get()));

        A.resize(coeffs.rows(), this->m_n_cols);

        //   b       eta                                         gamma
        A << coeffs, Vector::Constant(coeffs.rows(), eta_coeff), Vector::Zero(coeffs.rows());
        b = Vector::Constant(A.rows(), lower_bound);
    } else {
        Matrix tf = set.transformationMatrix(this->m_p);
        Matrix tf_FGamma = BRY_MATMUL(set.transformationMatrix(this->m_p), m_F_expec_Gamma_minus_I);
        Matrix coeffs = -BRY_MATMUL(this->getPhip(set.bernstein_deg_incr), tf_FGamma);
        //Matrix coeffs = -this->getPhip(set.bernstein_deg_incr) * set.transformationMatrix(this->m_p) * (m_F_expec_Gamma_minus_I);

        A.resize(coeffs.rows(), this->m_n_cols);

        //   b       eta                          gamma
        A << coeffs, Vector::Zero(coeffs.rows()), Vector::Ones(coeffs.rows());
        b = Vector::Zero(A.rows());
    }
    return std::make_pair(std::move(A), std::move(b));
}

#pragma once

#include "Problem.h"

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
const BRY::Matrix& BRY::PolyDynamicsProblem<DIM>::getPhim(bry_int_t bernstein_deg_incr) {
    auto it = m_Phi_m.find(bernstein_deg_incr);
    if (it == m_Phi_m.end()) {
        it = m_Phi_m.emplace(bernstein_deg_incr, BernsteinBasisTransform<DIM>::pwrToBernMatrix(this->barrier_deg, bernstein_deg_incr)).first;
    } 
    return it->second;
}

template <std::size_t DIM>
const BRY::Matrix& BRY::PolyDynamicsProblem<DIM>::getPhip(bry_int_t bernstein_deg_incr) {
    auto it = m_Phi_p.find(bernstein_deg_incr);
    if (it == m_Phi_p.end()) {
        it = m_Phi_p.emplace(bernstein_deg_incr, BernsteinBasisTransform<DIM>::pwrToBernMatrix(m_p, bernstein_deg_incr)).first;
    } 
    return it->second;
}

template <std::size_t DIM>
std::pair<BRY::Matrix, BRY::Vector> BRY::PolyDynamicsProblem<DIM>::calculateSetConstraints(ConstraintType constraint_type, const HyperRectangle<DIM>& set) {
    Matrix A;
    Vector b;
    bry_float_t lower_bound = 0.0;

    if (constraint_type != ConstraintType::Safe) {
        // Eta coeffs are 1 if the set is an initial set, otherwise they are zero
        bry_float_t eta_coeff = static_cast<bry_float_t>(constraint_type == ConstraintType::Init);
        // Lower bound is zero unless the set is of type unsafe, in which case the lb is 1
        lower_bound = static_cast<bry_float_t>(constraint_type == ConstraintType::Unsafe);
        // If the set is an initial set, then negate the b coefficients
        bry_float_t coeff_multiplier = (constraint_type == ConstraintType::Init) ? -1.0 : 1.0;

        Matrix coeffs = coeff_multiplier * this->getPhim(set.bernstein_deg_incr) * set.transformationMatrix(this->barrier_deg, filter.get());

        A.resize(coeffs.rows(), this->m_n_cols);

        //   b       eta                                         gamma
        A << coeffs, Vector::Constant(coeffs.rows(), eta_coeff), Vector::Zero(coeffs.rows());
        b = Vector::Constant(A.rows(), lower_bound);
    } else {
        Matrix coeffs = -this->getPhip(set.bernstein_deg_incr) * set.transformationMatrix(this->m_p) * (m_F_expec_Gamma_minus_I);

        A.resize(coeffs.rows(), this->m_n_cols);

        //   b       eta                          gamma
        A << coeffs, Vector::Zero(coeffs.rows()), Vector::Ones(coeffs.rows());
        b = Vector::Zero(A.rows());
    }
    return std::make_pair(std::move(A), std::move(b));
}

#pragma once

#include "LPSolver.h"
#include "MonomialFilter.h"
#include "Problem.h"

#include <list>
#include <memory>

namespace BRY {


template <std::size_t DIM>
struct SynthesisResult : public LPSolver::Result {
    SynthesisResult() = default;
    SynthesisResult(bry_int_t barrier_deg_) : barrier_deg(barrier_deg_) {}
    SynthesisResult(bry_int_t barrier_deg_, const std::shared_ptr<MonomialFilter<DIM>>& filter_) : barrier_deg(barrier_deg_), filter(filter_) {}

    BRY_INL bool isFilterApplied() const {return (bool)filter;}

    /// @brief Remove the current filter if one was supplied (otherwise do nothing)
    void removeFilter();

    bry_int_t barrier_deg = 0;

    std::shared_ptr<MonomialFilter<DIM>> filter = nullptr;
};

template <std::size_t DIM>
SynthesisResult<DIM> synthesize(const PolyDynamicsProblem<DIM>& problem, const std::string& solver_id = "clp");

template <std::size_t DIM>
SynthesisResult<DIM> synthesize(const ConstraintMatrices<DIM>& constraints, bry_int_t time_horizon, const std::string& solver_id = "clp");

template <std::size_t DIM>
SynthesisResult<DIM> synthesizeAdaptive(PolyDynamicsProblem<DIM> problem, bry_int_t max_iter, bry_int_t subdiv_per_iter, const std::string& solver_id = "clp");

void writeMatrixToFile(const Matrix& matrix, const std::string& filename);

}

#include "impl/Synthesis_impl.hpp"
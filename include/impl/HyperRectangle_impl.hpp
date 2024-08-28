#pragma once

#include "HyperRectangle.h"

#include "berry/Options.h"
#include "berry/MultiIndex.h"

#include <Eigen/Core>

#define COMP_DIFF_TOL 1e-3
//#define COMP_DIFF_TOL BRY_FLOAT_DIFF_TOL

template <std::size_t DIM>
BRY::HyperRectangle<DIM>::HyperRectangle(bry_float_t lower_default, bry_float_t upper_default)
    : lower_bounds(Eigen::Vector<bry_float_t, DIM>::Constant(lower_default))
    , upper_bounds(Eigen::Vector<bry_float_t, DIM>::Constant(upper_default))
{}

template <std::size_t DIM>
BRY::HyperRectangle<DIM>::HyperRectangle(bry_float_t lower_default, bry_float_t upper_default, bry_int_t degree_increase)
    : lower_bounds(Eigen::Vector<bry_float_t, DIM>::Constant(lower_default))
    , upper_bounds(Eigen::Vector<bry_float_t, DIM>::Constant(upper_default))
    , bernstein_deg_incr(degree_increase)
{}

template <std::size_t DIM>
bool BRY::HyperRectangle<DIM>::operator==(const HyperRectangle& other) const {
    return bernstein_deg_incr == other.bernstein_deg_incr
        && ((lower_bounds - other.lower_bounds).norm() < COMP_DIFF_TOL)
        && ((upper_bounds - other.upper_bounds).norm() < COMP_DIFF_TOL);
}

template <std::size_t DIM>
bool BRY::HyperRectangle<DIM>::operator<(const HyperRectangle& other) const {
    // Lexicographical comparison of two vectors
    auto lexLess = [](const auto& v1, const auto& v2) -> bool {
        for (bry_int_t i = 0; i < v1.size(); ++i) {
            if (v1[i] < (v2[i] - COMP_DIFF_TOL)) {
                return true;
            } else if (v1[i] > (v2[i] + COMP_DIFF_TOL)) {
                return false;
            }
        }
        return false;
    };

    if (bernstein_deg_incr < other.bernstein_deg_incr) {
        return true;
    } else if (bernstein_deg_incr > other.bernstein_deg_incr) {
        return false;
    } else if (lexLess(lower_bounds, other.lower_bounds)) {
        return true;
    } else if (lexLess(other.lower_bounds, lower_bounds)) {
        return false;
    } else if (lexLess(upper_bounds, other.upper_bounds)) {
        return true;
    } 
    return false;
}

template <std::size_t DIM>
BRY::Matrix BRY::HyperRectangle<DIM>::transformationMatrix(bry_int_t m, const MonomialFilter<DIM>* filter) const {

    Matrix T;
    if (!!filter) {
        T.resize(filter->nRemainingMonoms(), filter->nRemainingMonoms());
    } else {
        BRY::bry_int_t m_monoms = pow(m + 1, DIM);
        T.resize(m_monoms, m_monoms);
    }
    T.setZero();

    Eigen::Vector<bry_float_t, DIM> scale = scaleFromUnit();
    Eigen::Vector<bry_float_t, DIM> translation = translationFromUnit();

    const std::vector<bool>* filter_flags = !filter ? nullptr : &filter->flags();

    for (auto row_midx = mIdxW(DIM, m + 1); !row_midx.last(); ++row_midx) {
        if (!!filter && (*filter_flags)[row_midx.inc().wrappedIdx()]) {
            continue;
        }

        std::vector<bry_int_t> index_bounds(row_midx.size());
        for (std::size_t d = 0; d < DIM; ++d) {
            index_bounds[d] = m + 1 - row_midx[d];
        }

        for (auto col_midx = mIdxBEW(index_bounds, m + 1); !col_midx.last(); ++col_midx) {
            if (!!filter && (*filter_flags)[row_midx.inc().wrappedIdx() + col_midx.inc().wrappedIdx()]) {
                continue;
            }

            bry_float_t element = 1.0;
            for (std::size_t j = 0; j < DIM; ++j) {
                element *= binom(row_midx[j] + col_midx[j], row_midx[j]);
                element *= std::pow(scale[j], row_midx[j]);
                element *= std::pow(translation[j], col_midx[j]);
            }

            bry_int_t r = row_midx.inc().wrappedIdx();
            bry_int_t c = row_midx.inc().wrappedIdx() + col_midx.inc().wrappedIdx();
            if (!!filter) {
                r = filter->newWrappedIdx(r);
                c = filter->newWrappedIdx(c);
            }

            T(r, c) = element;
        }
    }
    return T;
}

template <std::size_t DIM>
Eigen::Vector<BRY::bry_float_t, DIM> BRY::HyperRectangle<DIM>::scaleToUnit() const {
    Eigen::Vector<BRY::bry_float_t, DIM> diff = upper_bounds - lower_bounds;
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT((diff.array() > 0.0).all(), "Upper bounds are not greater than lower bounds");
#endif
    return diff.array().pow(-1.0);
}

template <std::size_t DIM>
Eigen::Vector<BRY::bry_float_t, DIM> BRY::HyperRectangle<DIM>::translationToUnit() const {
    return -lower_bounds.cwiseProduct(scaleToUnit());
}

template <std::size_t DIM>
Eigen::Vector<BRY::bry_float_t, DIM> BRY::HyperRectangle<DIM>::scaleFromUnit() const {
    Eigen::Vector<BRY::bry_float_t, DIM> diff = upper_bounds - lower_bounds;
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT((diff.array() > 0.0).all(), "Upper bounds are not greater than lower bounds");
#endif
    return diff;
}

template <std::size_t DIM>
Eigen::Vector<BRY::bry_float_t, DIM> BRY::HyperRectangle<DIM>::translationFromUnit() const {
    return lower_bounds;
}

template <std::size_t DIM>
std::vector<BRY::HyperRectangle<DIM>> BRY::HyperRectangle<DIM>::subdivide(uint32_t integer_division) const {
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT(integer_division > 0, "Integer division must be greater than 0");
#endif

    std::vector<BRY::HyperRectangle<DIM>> subsets;
    subsets.reserve(pow(integer_division, DIM));
    Eigen::Vector<bry_float_t, DIM> diff = (upper_bounds - lower_bounds) / static_cast<float>(integer_division);
    for (auto midx = mIdx(DIM, integer_division); !midx.last(); ++midx) {
        Eigen::Vector<bry_float_t, DIM> shift = Eigen::Vector<bry_float_t, DIM>::Zero();
        for (std::size_t d = 0; d < DIM; ++d) {
            shift[d] = midx[d] * diff[d];
        }
        BRY::HyperRectangle<DIM> subset;
        subset.lower_bounds = lower_bounds + shift;
        subset.upper_bounds = subset.lower_bounds + diff;
        subsets.push_back(std::move(subset));
    }
    return subsets;
}

//template <std::size_t DIM>
//std::vector<BRY::HyperRectangle<DIM>> BRY::HyperRectangle<DIM>::subdivideByPercent(const Eigen::Vector<bry_float_t, DIM>& division_percent_point) const {
//
//#ifdef BRY_ENABLE_BOUNDS_CHECK
//    for (auto pc : division_percent_point) {
//        ASSERT(percent_point > 0.0 && percent_point < 1.0, "Split volume percentage must be between 0.0 and 1.0");
//    }
//#endif
//
//    std::vector<BRY::HyperRectangle<DIM>> subsets;
//    subsets.reserve(pow(integer_division, DIM));
//    Eigen::Vector<bry_float_t, DIM> to_division_point = (upper_bounds - lower_bounds).array() * 
//    for (auto midx = mIdx(DIM, 2); !midx.last(); ++midx) {
//        Eigen::Vector<bry_float_t, DIM> shift = Eigen::Vector<bry_float_t, DIM>::Zero();
//        for (std::size_t d = 0; d < DIM; ++d) {
//            shift[d] = midx[d] * diff[d];
//        }
//        BRY::HyperRectangle<DIM> subset;
//        subset.lower_bounds = lower_bounds + shift;
//        subset.upper_bounds = subset.lower_bounds + diff;
//        subsets.push_back(std::move(subset));
//    }
//    return subsets;
//}

template <std::size_t DIM>
std::pair<BRY::HyperRectangle<DIM>, BRY::HyperRectangle<DIM>> BRY::HyperRectangle<DIM>::split(bry_int_t dim) const {
    return splitByPercent(dim, 0.5);
}

template <std::size_t DIM>
std::pair<BRY::HyperRectangle<DIM>, BRY::HyperRectangle<DIM>> BRY::HyperRectangle<DIM>::split(bry_int_t dim, bry_float_t split_point) const {
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT(split_point > lower_bounds[dim] && split_point < upper_bounds[dim], "Split point is not contained within hyperrectangle along split dimension");
#endif
    std::pair<HyperRectangle<DIM>, HyperRectangle<DIM>> split_pieces = std::make_pair(*this, *this);
    split_pieces.first.upper_bounds[dim] = split_point;
    split_pieces.second.lower_bounds[dim] = split_point;
    return split_pieces;
}

template <std::size_t DIM>
std::pair<BRY::HyperRectangle<DIM>, BRY::HyperRectangle<DIM>> BRY::HyperRectangle<DIM>::splitByPercent(bry_int_t dim, bry_float_t percent_point) const {
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT(percent_point > 0.0 && percent_point < 1.0, "Split volume percentage must be between 0.0 and 1.0");
#endif
    std::pair<HyperRectangle<DIM>, HyperRectangle<DIM>> split_pieces = std::make_pair(*this, *this);

    bry_float_t split_point = lower_bounds[dim] + percent_point * (upper_bounds[dim] - lower_bounds[dim]);
    split_pieces.first.upper_bounds[dim] = split_point;
    split_pieces.second.lower_bounds[dim] = split_point;
    return split_pieces;

}
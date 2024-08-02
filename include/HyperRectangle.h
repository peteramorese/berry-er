#pragma once

#include "berry/Options.h"

namespace BRY {

template <std::size_t DIM>
struct HyperRectangle {
    public:
        Eigen::Vector<bry_float_t, DIM> lower_bounds = Eigen::Vector<bry_float_t, DIM>::Constant(0.0);
        Eigen::Vector<bry_float_t, DIM> upper_bounds = Eigen::Vector<bry_float_t, DIM>::Constant(1.0);

    public:
        HyperRectangle() = default;

        /// @brief Construct with default uppper and lower bounds (same for all dimensions)
        /// @param lower_default Lower bound for all dims
        /// @param upper_default Upper bound for all dims
        HyperRectangle(bry_float_t lower_default, bry_float_t upper_default);

        Matrix transformationMatrix(bry_int_t m) const;

        BRY_INL Eigen::Vector<bry_float_t, DIM> scaleToUnit() const;
        BRY_INL Eigen::Vector<bry_float_t, DIM> translationToUnit() const;
        BRY_INL Eigen::Vector<bry_float_t, DIM> scaleFromUnit() const;
        BRY_INL Eigen::Vector<bry_float_t, DIM> translationFromUnit() const;

        std::vector<HyperRectangle<DIM>> subdivide(uint32_t integer_division = 2) const;

        std::vector<HyperRectangle<DIM>> subdivideByPercent(const Eigen::Vector<bry_float_t, DIM>& division_percent_point) const;

        /// @brief Split in half
        /// @param dim Dimension to split along
        /// @return Two halves of original hyperrectangle
        std::pair<HyperRectangle<DIM>, HyperRectangle<DIM>> split(bry_int_t dim) const;

        /// @brief Unevenly split in two
        /// @param dim Dimension to split along
        /// @param split_point Coordinate along splitting dimension to divide the rectangles at
        /// @return Split hyperrectangles
        std::pair<HyperRectangle<DIM>, HyperRectangle<DIM>> split(bry_int_t dim, bry_float_t split_point) const;

        /// @brief Unevenly split in two by percent volume
        /// @param dim Dimension to split along
        /// @param percent Number in `(0.0, 1.0)`. `percent` volume in first partition, `1-percent` voluming in second
        /// @return Split hyperrectangles
        std::pair<HyperRectangle<DIM>, HyperRectangle<DIM>> splitByPercent(bry_int_t dim, bry_float_t percent_point) const;
};

}

#include "impl/HyperRectangle_impl.hpp"
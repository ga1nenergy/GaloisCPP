#ifndef MISC_HPP
#define MISC_HPP

#include <type_traits>
#include <algorithm>
#include "galoiscpp.h"

template <class T>
bool is_degraded(const std::vector<std::vector<T>> &matrix ) {
    std::vector<T> zero_row(matrix[0].size());
    if (std::is_same<T, galoiscpp::GFelement>::value) {
        std::fill(zero_row.begin(), zero_row.end(), galoiscpp::GFelement(matrix[0][0].getField(), 0));
    } else {
        std::fill(zero_row.begin(), zero_row.end(), 0);
    }

    return std::find(matrix.begin(), matrix.end(), zero_row) != matrix.end();
}

template <class T>
std::vector<std::vector<T>> drop_zero_rows(const std::vector<std::vector<T>> &matrix) {
    auto new_matrix = matrix;
    std::vector<T> zero_row(matrix[0].size());

    if (std::is_same<T, galoiscpp::GFelement>::value) {
        std::fill(zero_row.begin(), zero_row.end(), galoiscpp::GFelement(matrix[0][0].getField(), 0));
    } else {
        std::fill(zero_row.begin(), zero_row.end(), 0);
    }

    for (auto i = std::find(new_matrix.begin(), new_matrix.end(), zero_row); i != new_matrix.end();
                i = std::find(new_matrix.begin(), new_matrix.end(), zero_row)) {
        auto dist = std::distance(new_matrix.begin(), i);
//        std::cout << "dist: " << dist << std::endl;
        new_matrix.erase(i);
        std::for_each(new_matrix.begin(), new_matrix.end(), [dist](auto &row) {
           row.erase(row.begin() + dist);
        });
    }

    return new_matrix;
}

#endif
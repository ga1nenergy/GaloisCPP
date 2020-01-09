//
// Created by ga1nenergy on 23.07.2019.
//

#include "gfalgorithms.h"

#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>

namespace galoiscpp {
    /* TODO
     * 1) add vector elements' fields compatibility check
     */
    GFpoly newton_identitites(const std::vector<GFelement> &pow_sum) {
        const GaloisField *field = pow_sum[0].getField();
        size_t t = pow_sum.size();
        std::vector<GFelement> sigma = {GFelement(field, 1)};

        for (int i = 1; i <= t; i++) {
            std::vector<GFelement> rev(pow_sum.begin(), pow_sum.begin() + i);
            std::reverse(rev.begin(), rev.end());

            auto left_side = GFelement::dot(sigma, rev);
            auto res = -left_side.summed_times(i);
            sigma.push_back(res);

//            for (int j = 0; j < field->get_size(); j++) {
//                GFelement cur_sigma(field, j);
//                auto buf = -cur_sigma;
//                if (left_side == -cur_sigma.sum_times(i)) {
//                    sigma.push_back(cur_sigma);
//                    break;
//                }
//            }
        }

        while (sigma[sigma.size() - 1] == 0) {
            sigma.pop_back();
        }

        GFpoly res(sigma);
        return res;
    }

    /* TODO
     * 1) improve logic
     */
    std::vector<std::vector<Fint>> gauss_method_modulus(const std::vector<std::vector<Fint>> &matrix, Fint modulus) {
        // Forward
        auto transformed_matrix = matrix;
        size_t r = transformed_matrix.size();
        size_t n = transformed_matrix[0].size();
        size_t k = n - r;

        std::vector<Fint> zero_row(n, 0);
        auto padding = 0;

        for (size_t i = 0; (i + padding) < transformed_matrix.size(); i++) {
//            auto coef = field.inverse(transformed_matrix[i][i]);
            Fint coef;
            if (transformed_matrix[i][i + padding] == 0) {
                auto flag = 1;
                for (auto j = i + 1; j < transformed_matrix.size(); j++) {
                    if (transformed_matrix[j][i] != 0) {
                        transformed_matrix[j].swap(transformed_matrix[i]);
                        flag = 0;
                        break;
                    }
                }

                if (flag) {
                    padding++;
                    continue;
                }
            }
            for (size_t z = 1; z < modulus; z++) {
                auto buf = (z * transformed_matrix[i][i]) % (modulus);
                if ((z * transformed_matrix[i][i]) % (modulus) == 1) {
                    coef = z;
                    break;
                }
            }

            for (size_t j = i + padding; j < n; j++) {
//                transformed_matrix[i][j] = field.multiply(coef, transformed_matrix[i][j]);
                if (transformed_matrix[i][j] != 0) {
                    transformed_matrix[i][j] = (coef * transformed_matrix[i][j]) % (modulus);
                }
            }

//            std::cout << std::endl;
//            std::copy(transformed_matrix[i].begin(), transformed_matrix[i].end(), std::ostream_iterator<Fint>(std::cout, " "));
//            std::cout << std::endl;

            for (size_t j = i + 1 + padding; j < transformed_matrix.size(); j++) {
                coef = transformed_matrix[j][i];
//                std::cout << "Row " << j << std::endl;
                if (transformed_matrix[j][i] != 0) {
                    for (size_t p = i; p < n; p++) {
                        //                    transformed_matrix[j][p] = field.subtract(transformed_matrix[j][p], field.multiply(coef,
                        //                                                                                        transformed_matrix[i][p]));
                        auto mul = 0;
                        if (transformed_matrix[i][p] != 0)
                            mul = (coef * transformed_matrix[i][p]) % (modulus);
//                                            std::cout << coef << "*" << transformed_matrix[i][p] << "=" << mul << std::endl;
                        auto buf = transformed_matrix[j][p] - mul;
                        if (buf < 0) buf = modulus + buf;
//                                            std::cout << transformed_matrix[j][p] << "-" << mul << "=" << buf << std::endl;

                        transformed_matrix[j][p] = buf;
                    }
                }
//                std::copy(transformed_matrix[j].begin(), transformed_matrix[j].end(), std::ostream_iterator<Fint>(std::cout, " "));
//                std::cout << std::endl;
            }

//            for (auto f = std::find(transformed_matrix.begin(), transformed_matrix.end(), zero_row); f != transformed_matrix.end();
//                 f = std::find(transformed_matrix.begin(), transformed_matrix.end(), zero_row)) {
//                transformed_matrix.erase(f);
//            }

//            for (auto & row : transformed_matrix) {
//                for (auto & e : row) {
//                    std::cout << e << " ";
//                }
//                std::cout << std::endl;
//            }
//            std::cout << std::endl;
        }
//        std::cout << std::endl;

        // Backward
        for (int i = transformed_matrix.size() - 1; i >= 1; i--) {
            for (int j = i - 1; j >= 0; j--){
                auto coef = transformed_matrix[j][i];
                if (transformed_matrix[j][i] != 0) {
                    for (int p = i; p < n; p++) {
//                        transformed_matrix[j][p] = field.subtract(transformed_matrix[j][p], field.multiply(coef,
//                                                                                            transformed_matrix[i][p]));
                        auto mul = 0;
                        if (transformed_matrix[i][p] != 0)
                            mul = (coef * transformed_matrix[i][p]) % (modulus);
                        transformed_matrix[j][p] = transformed_matrix[j][p] - mul;
                        if (transformed_matrix[j][p] < 0) transformed_matrix[j][p] = modulus + transformed_matrix[j][p];
                    }
                }
            }

//            for (auto & row : transformed_matrix) {
//                for (auto & e : row) {
//                    std::cout << e << " ";
//                }
//                std::cout << std::endl << std::endl;
//            }
        }

//        for (auto & row : transformed_matrix) {
//            for (auto & e : row) {
//                std::cout << e << " ";
//            }
//            std::cout << std::endl;
//        }

        return transformed_matrix;
    }

    /* TODO
 * 1) improve logic
 */
    std::vector<std::vector<GFelement>> gauss_method(const std::vector<std::vector<GFelement>> &matrix) {
        bool DEBUG = false;

        std::vector<int> permutations(matrix.size());
        std::iota(permutations.begin(), permutations.end(), 0);

        size_t r = matrix.size();
        size_t n = matrix[0].size();

        auto transformed_matrix = matrix;

        // Forward
        for (size_t i = 0; i < r; i++) {
            if (transformed_matrix[i][i] == 0) {
                for (auto j = i + 1; j < transformed_matrix.size(); j++) {
                    if (transformed_matrix[j][i] != 0) {
                        iter_swap(transformed_matrix.begin() + i, transformed_matrix.begin() + j);
                        iter_swap(permutations.begin() + i, permutations.begin() + j);
                    }
                }
            }

            auto coef = transformed_matrix[i][i].inverse();
            if (coef != 0) {
                for (size_t j = i; j < n; j++) {
                    transformed_matrix[i][j] = transformed_matrix[i][j] * coef;
                }
            }

//            for (size_t j = i + 1; j < r; j++) {
//                coef = transformed_matrix[j][i];
//                for (size_t p = i; p < n; p++) {
//                    if (coef != 0) {
//                        transformed_matrix[j][p] = transformed_matrix[j][p] - coef * transformed_matrix[i][p];
//                    }
//                }
//            }

            for (size_t j = i + 1; j < r; j++) {
                coef = transformed_matrix[j][i];
                for (size_t p = i; p < n; p++) {
                    if (transformed_matrix[i][p] != 0)

                    if (coef != 0) {
                        transformed_matrix[j][p] = transformed_matrix[j][p] - coef * transformed_matrix[i][p];
                    }
                }
            }

            if (DEBUG) {
                std::cout << std::endl << "Fwd. Step " << i << ": " << std::endl;
                for (auto &row : transformed_matrix) {
                    for (auto &e : row) {
                        std::cout << e << " ";
                    }
                    std::cout << std::endl;
                }
            }
        }

        if (DEBUG) {
            std::cout << std::endl << "Forward: " << std::endl;
            for (auto &row : transformed_matrix) {
                for (auto &e : row) {
                    std::cout << e << " ";
                }
                std::cout << std::endl;
            }
        }

        // Backward
        for (int i = r - 1; i >= 0; i--) {
            for (int j = i - 1; j >= 0; j--){
                auto coef = transformed_matrix[j][i];
                if (coef != 0) {
                    for (int p = i; p < n; p++) {
                        transformed_matrix[j][p] = transformed_matrix[j][p] - coef * transformed_matrix[i][p];
                    }
                }
            }

            if (DEBUG) {
                std::cout << std::endl << "Bkw. Step " << i << ": " << std::endl;
                for (auto &row : transformed_matrix) {
                    for (auto &e : row) {
                        std::cout << e << " ";
                    }
                    std::cout << std::endl;
                }
            }
        }

        if (DEBUG) {
            std::cout << std::endl << "Backward: " << std::endl;
            for (auto &row : transformed_matrix) {
                for (auto &e : row) {
                    std::cout << e << " ";
                }
                std::cout << std::endl;
            }
        }

        return transformed_matrix;
    }

    std::vector<std::vector<GFelement>> gauss_method_mixed(const std::vector<std::vector<GFelement>> &matrix) {
        size_t r = matrix.size();
        size_t n = matrix[0].size();

        auto transformed_matrix = matrix;

        // Forward
        for (size_t i = 0; i < r; i++) {
            auto coef = transformed_matrix[i][i].inverse();
            for (size_t j = i; j < n; j++) {
                transformed_matrix[i][j] = transformed_matrix[i][j] * coef;
            }

            for (size_t j = i + 1; j < r; j++) {
                coef = transformed_matrix[j][i];
                for (size_t p = i; p < n; p++) {
                    transformed_matrix[j][p] = transformed_matrix[j][p] - coef * transformed_matrix[i][p];
                }
            }
        }

//        for (auto & row : transformed_matrix) {
//            for (auto &e : row) {
//                std::cout << e << " ";
//            }
//            std::cout << std::endl;
//        }

        // Backward
        for (int i = r - 1; i >= 0; i--) {
            for (int j = i - 1; j >= 0; j--){
                auto coef = transformed_matrix[j][i];
                for (int p = i; p < n; p++) {
                    transformed_matrix[j][p] = transformed_matrix[j][p] - coef * transformed_matrix[i][p];
                }
            }
        }

//        for (auto & row : transformed_matrix) {
//            for (auto &e : row) {
//                std::cout << e << " ";
//            }
//            std::cout << std::endl;
//        }

        return transformed_matrix;
    }

    GFpoly lagrange_interpolation(const std::vector<GFelement> &x, const std::vector<GFelement> &y) {
        auto DEBUG = false;
        const GaloisField *gf = y[0].getField();

        GFpoly f_x(gf, 0); f_x[0] = 0;
        for (auto i = 0; i < x.size(); i++) {
            if (DEBUG) {
                std::cout << "Start of iter #" << i << std::endl;
                std::cout << "f_x: " << f_x << std::endl;
            }

            GFpoly prod(gf, 0); prod[0] = 1;
            for (auto j = 0; j < y.size(); j++) {
                if (x[i] != x[j]) {
                    GFpoly buf(gf, 1);
                    buf[0] = -x[j]; buf[1] = 1;
                    prod = prod * (buf / (x[i] - x[j]));
                }

                if (DEBUG) {
                    std::cout << "Prod: " << prod << std::endl;
                }
            }
            f_x = f_x + y[i] * prod;
            if (DEBUG) {
                std::cout << "End of iter #" << i << std::endl;
                std::cout << "f_x: " << f_x << std::endl;
            }
        }

        return f_x;
    }
}
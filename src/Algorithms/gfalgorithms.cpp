//
// Created by ga1nenergy on 23.07.2019.
//

#include "gfalgorithms.h"

#include <vector>
#include <algorithm>
#include <iostream>

namespace galoiscpp {
    /* TODO
     * 1) add vector elements' fields compatibility check
     */
    GFpoly newton_identitites(const std::vector<GFelement> &pow_sum) {
        GaloisField field = pow_sum[0].getField();
        size_t t = pow_sum.size();
        std::vector<GFelement> sigma = {GFelement(field, 1)};

        for (int i = 1; i <= t; i++) {
            std::vector<GFelement> rev(pow_sum.begin(), pow_sum.begin() + i);
            std::reverse(rev.begin(), rev.end());

            auto left_side = GFelement::dot(sigma, rev);
            auto res = (-left_side).summed_times(i);
            sigma.push_back(res);

//            for (int j = 0; j < field.get_size(); j++) {
//                GFelement cur_sigma(field, j);
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
        size_t r = matrix.size();
        size_t n = matrix[0].size();
        size_t k = n - r;

        auto transformed_matrix = matrix;

        for (size_t i = 0; i < r; i++) {
//            auto coef = field.inverse(transformed_matrix[i][i]);
            Fint coef;
            for (size_t z = 1; z < modulus; z++) {
                auto buf = (z * transformed_matrix[i][i]) % (modulus);
                if ((z * transformed_matrix[i][i]) % (modulus) == 1) {
                    coef = z;
                    break;
                }
            }

            for (size_t j = i; j < n; j++) {
//                transformed_matrix[i][j] = field.multiply(coef, transformed_matrix[i][j]);
                if (transformed_matrix[i][j] != 0) {
                    transformed_matrix[i][j] = (coef * transformed_matrix[i][j]) % (modulus);
                }
            }

            for (size_t j = i + 1; j < r; j++) {
                coef = transformed_matrix[j][i];
//                std::cout << "Row " << j << std::endl;
                if (transformed_matrix[j][i] != 0) {
                    for (size_t p = i; p < n; p++) {
                        //                    transformed_matrix[j][p] = field.subtract(transformed_matrix[j][p], field.multiply(coef,
                        //                                                                                        transformed_matrix[i][p]));
                        auto mul = 0;
                        if (transformed_matrix[i][p] != 0)
                            mul = (coef * transformed_matrix[i][p]) % (modulus);
                        //                    std::cout << coef << "*" << transformed_matrix[i][p] << "=" << mul << std::endl;
                        auto buf = transformed_matrix[j][p] - mul;
                        if (buf < 0) buf = modulus + buf;
                        //                    std::cout << transformed_matrix[j][p] << "-" << mul << "=" << buf << std::endl;

                        transformed_matrix[j][p] = buf;
                    }
                }
            }

//            for (auto & row : transformed_matrix) {
//                for (auto & e : row) {
//                    std::cout << e << " ";
//                }
//                std::cout << std::endl;
//            }
//            std::cout << std::endl;
        }

        // Backward
        for (int i = r - 1; i >= 1; i--) {
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
                if (coef != 0) {
                    for (size_t p = i; p < n; p++) {
                        transformed_matrix[j][p] = transformed_matrix[j][p] - coef * transformed_matrix[i][p];
                    }
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
                if (coef != 0) {
                    for (int p = i; p < n; p++) {
                        transformed_matrix[j][p] = transformed_matrix[j][p] - coef * transformed_matrix[i][p];
                    }
                }
            }
        }

//        for (auto & row : transformed_matrix) {
//            for (auto & e : row) {
//                std::cout << e << " ";
//            }
//            std::cout << std::endl;
//        }

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
}
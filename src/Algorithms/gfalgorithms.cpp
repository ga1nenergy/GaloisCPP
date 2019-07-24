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
        GaloisField *field = pow_sum[1].getField();
        size_t t = pow_sum.size();
        std::vector<GFelement> sigma = {GFelement(field, 1)};

        for (int i = 1; i <= t; i++) {
            std::vector<GFelement> rev(pow_sum.begin(), pow_sum.begin() + i);
            std::reverse(rev.begin(), rev.end());

            auto left_side = GFelement::dot(sigma, rev);

            for (int j = 0; j < field->get_size(); j++) {
                GFelement cur_sigma(field, j);
                if (left_side == -cur_sigma.sum_times(i)) {
                    sigma.push_back(cur_sigma);
                    break;
                }
            }
        }

        while (sigma[sigma.size() - 1] == 0) {
            sigma.pop_back();
        }

        GFpoly res(sigma);
        return res;
    }

    std::vector<std::vector<Fint>> gauss_method_modulus(const std::vector<std::vector<Fint>> &matrix, Fint modulus) {
        // Forward 
    }
}
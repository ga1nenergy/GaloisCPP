//
// Created by ga1nenergy on 23.07.2019.
//

#include <vector>
#include "gfelement.h"
#include "gfpoly.h"

#ifndef GALOISCPP_REMAINDER_CODE_H
#define GALOISCPP_REMAINDER_CODE_H

using namespace galoiscpp;

namespace remainder_code {
    struct scenario_t {
        int n_all, n_strong;
        int t_weak, t_strong;
    };

    std::vector<std::vector<GFelement>> parity_check_matrix(const std::vector<GFelement> &locators, int t);
    std::vector<std::vector<GFelement>> parity_check_matrix_punct(const std::vector<GFelement> &locators, scenario_t sc);

    GFpoly find_remainder(const std::vector<int> &msg, const std::vector<std::vector<GFelement>> &H);
    std::vector<GFpoly> find_strong_remainders(const std::vector<int> &msg, const std::vector<GFelement> &locators_strong,
                                               const std::vector<std::vector<GFelement>> &H,
                                               const std::vector<std::vector<GFelement>> &combs);
    std::vector<std::vector<GFelement>> find_locator_combinations(const std::vector<GFelement> &locators, int len);
    void locator_combination_iter(const std::vector<GFelement> &locators, std::vector<std::vector<GFelement>> &combs,
                                  std::vector<GFelement> &counter, int elem_number);

    std::vector<int> decode(const std::vector<int> &received_msg, const std::vector<GFelement> &locators,
                            const std::vector<std::vector<GFelement>> &H, const GFpoly &rem);
    std::vector<std::vector<int>> decode_strong(const std::vector<int> &received_msg, const std::vector<GFelement> &locators,
                                                const std::vector<std::vector<GFelement>> &H, const GFpoly &rem,
                                                const std::vector<std::vector<GFelement>> &combs);
    std::vector<int> decode_multiuser(const std::vector<int> &received_msg, const std::vector<GFelement> &locators,
                                      const std::vector<std::vector<GFelement>> &H, const GFpoly &rem,
                                      const std::vector<std::vector<GFelement>> &combs, scenario_t sc);

    std::vector<int> decode_multiuser_separable(const std::vector<int> &received_msg, const std::vector<GFelement> &locators,
                                                const std::vector<std::vector<GFelement>> &H,
                                                const std::vector<GFelement> &syndrome, scenario_t sc);
};

#endif //GALOISCPP_REMAINDER_CODE_H

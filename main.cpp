#include <GFpoly/representations.h>
#include "gfelement.h"
#include "galoisfield.h"
#include "gfpoly.h"
#include "remainder_code.h"

#include <iostream>
#include <vector>
#include <algorithm>

using namespace galoiscpp;

int main() {
    int n = 12, n_sub = 6;
    int p = 7, m = 2;
    int t_strong = 4, t_weak = 2;

    remainder_code::scenario_t sc = {n, n_sub, t_weak, t_strong};

    GaloisField field(p, m);

    std::vector<Fint> subfield = field.find_subfield(1);

    std::vector<GFelement> subfield_locators = GFelement::to_gf(&field, subfield);
    subfield_locators.resize(n_sub);

    std::vector<GFelement> locators = subfield_locators;
    for (int i = 1; i < field.get_size() && locators.size() < n; i++) {
        if (std::find(subfield.begin(), subfield.end(), i) == subfield.end())
            locators.emplace_back(&field, i);
    }

    std::vector<int> msg = {1, 2, 2, 1, 2, 1, 0, 0, 0, 0, 0, 0};
//    msg.resize(sc.n_strong); //!!!
//    std::vector<std::vector<GFelement>> H = remainder_code::parity_check_matrix(locators, t);
    std::vector<std::vector<GFelement>> H = remainder_code::parity_check_matrix_punct(locators, sc);

    std::vector<GFelement> syndrome(H.size());
    for (size_t i = 0; i < syndrome.size(); i++) {
        syndrome[i] = GFelement::dotint(H[i], msg);
    }

//    for (auto & row : H) {
//        row.resize(sc.n_strong);
//    }

    for (auto & row : H) {
        for (auto & e : row) {
            std::cout << e << " ";
        }
        std::cout << std::endl;
    }

    GFpoly rem = remainder_code::find_remainder(msg, H);
    std::cout << "Remainder: " << rem << std::endl;

    std::vector<int> msg_error = {2, 2, 2, 1, 2, 1, 0, 0, 0, 0, 0, 0};  //! Debug this case. Error in root() function

//    msg_error.resize(sc.n_strong);  //!!!
//    auto decoded_msg = remainder_code::decode(msg_error, locators, H, rem);

//    auto combs = remainder_code::find_locator_combinations(subfield_locators, sc.t_weak);
//    auto decoded_msg = remainder_code::decode_multiuser(msg_error, locators, H, rem, combs, sc);
    auto decoded_msg = remainder_code::decode_multiuser_separable(msg_error, locators, H, syndrome, sc);

    std::cout << "-----------------------------------" << std::endl;
    std::cout << "Decoded_codeword: ";
    for (auto const & e : decoded_msg) {
        std::cout << e << " ";
    }
    std::cout << std::endl;

    return 0;
}
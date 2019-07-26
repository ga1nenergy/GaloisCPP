//#include <Codec/representations.h>
//#include "gfelement.h"
//#include "galoisfield.h"
//#include "gfpoly.h"
//#include "remainder_code.h"
//#include <Algorithms/gfalgorithms.h>

#include "api/galoiscpp.h"

#include <iostream>
#include <vector>
#include <algorithm>

using namespace galoiscpp;

int main() {
    int n = 12, n_strong = 6;
    int p = 7, m = 2;
    int t_strong = 4, t_weak = 2;

    GaloisField field(p, m); std::cout << field << std::endl;

    std::vector<Fint> subfield = field.find_subfield(1);

    std::vector<GFelement> subfield_locators = GFelement::to_gf(&field, subfield);
    subfield_locators.resize(n_strong);

    std::vector<GFelement> locators = subfield_locators;
    std::vector<GFelement> field_locators;
    for (int i = 1; i < field.get_size() && locators.size() < n; i++) {
        if (std::find(subfield.begin(), subfield.end(), i) == subfield.end()) {
            locators.emplace_back(&field, i);
        }
    }

    for (int i = 1; i < field.get_size(); i++) {
        if (std::find(subfield.begin(), subfield.end(), i) == subfield.end())
            field_locators.emplace_back(&field, i);
    }

    std::vector<Fint> msg = {1, 1, 2, 1, 0, 0, 1, 0, 1, 0, 0, 0};
    std::vector<Fint> msg_strong(msg.begin(), msg.begin() + n_strong);

    GFpoly g_x(&field, t_strong + 1); g_x[g_x.getDegree()] = 1;
    GFpoly g_x_weak(&field, t_weak + 1); g_x_weak[g_x_weak.getDegree()] = 1;
    auto field_locators_zero = field_locators; field_locators_zero.emplace_back(&field, 0);
    auto g_x_l_x = roots_to_poly(field_locators_zero);

    auto poly = locator_polynomial(subfield_locators, msg_strong);
    auto qr = poly / g_x_l_x;
    qr = qr[1] / g_x;
    GFpoly rem_strong = qr[1];
    std::cout << "Strong remainder: " << rem_strong << std::endl;

    poly = locator_polynomial(locators, msg);
    qr = poly / g_x_weak;
    GFpoly rem_weak = qr[1];
    std::cout << "Weak remainder: " << rem_weak << std::endl;

    auto rem_common = GFpoly::convmod(rem_strong, rem_weak, g_x);
    std::cout << "Common remainder: " << rem_common << std::endl;

    std::vector<Fint> msg_error = {1, 1, 2, 1, 0, 0, 1, 0, 1, 1, 0, 0};
    msg_strong = std::vector<Fint>(msg_error.begin(), msg_error.begin() + n_strong);
    poly = locator_polynomial(subfield_locators, msg_strong);
    qr = poly / g_x_l_x;
    std::cout << "Received strong remainder: " << qr[1] << std::endl;
    qr = qr[1] / g_x;
    GFpoly rem_strong_received = qr[1];
    std::cout << "Received strong remainder: " << rem_strong_received << std::endl;

    auto u_strong = GFpoly::deconvmod(rem_strong_received, rem_strong, g_x);
    std::cout << "Strong inc error poly: " << u_strong << std::endl;

    return 0;
}
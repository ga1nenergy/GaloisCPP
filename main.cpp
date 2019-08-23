#include "galoiscpp.h"

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

    std::vector<GFelement> subfield_locators = GFelement::to_gf(field, subfield);
    subfield_locators.resize(n_strong);

    std::vector<GFelement> locators = subfield_locators;
    std::vector<GFelement> field_locators;
    for (int i = 1; i < field.get_size() && locators.size() < n; i++) {
        if (std::find(subfield.begin(), subfield.end(), i) == subfield.end()) {
            locators.emplace_back(field, i);
        }
    }

    for (int i = 1; i < field.get_size(); i++) {
        if (std::find(subfield.begin(), subfield.end(), i) == subfield.end())
            field_locators.emplace_back(field, i);
    }

    std::vector<int> msg = {1, 1, 2, 1, 0, 0, 1, 0, 1, 0, 0, 0};
    std::vector<int> msg_strong(msg.begin(), msg.begin() + n_strong);

    GFpoly g_x(field, t_strong + 1); g_x[g_x.getDegree()] = 1;
    GFpoly g_x_weak(field, t_weak + 1); g_x_weak[g_x_weak.getDegree()] = 1;
    auto field_locators_zero = field_locators; field_locators_zero.emplace_back(field, 0);
    auto g_x_l_x = roots_to_poly(field_locators_zero);

    auto poly = locator_polynomial(subfield_locators, msg_strong);
    GFpoly rem_strong = GFpoly::deconvmod(poly, g_x_l_x, g_x);
    std::cout << "Strong remainder: " << rem_strong << std::endl;

    poly = locator_polynomial(locators, msg);
    GFpoly rem_weak = GFpoly::deconvr(poly, g_x_weak);
    std::cout << "Weak remainder: " << rem_weak << std::endl;

    auto rem_common = GFpoly::convmod(rem_strong, rem_weak, g_x);
    std::cout << "Common remainder: " << rem_common << std::endl;

    std::vector<int> msg_error = {1, 1, 2, 1, 0, 0, 1, 0, 1, 1, 0, 0};
    msg_strong = std::vector<int>(msg_error.begin(), msg_error.begin() + n_strong);
    poly = locator_polynomial(subfield_locators, msg_strong);
    std::cout << "Received strong remainder: " << GFpoly::deconvr(poly, g_x_l_x) << std::endl;
    GFpoly rem_strong_received = GFpoly::deconvmod(poly, g_x_l_x, g_x);
    std::cout << "Received strong remainder: " << rem_strong_received << std::endl;

    auto u_strong = GFpoly::deconvmod(rem_strong_received, rem_strong, g_x);
    std::cout << "Strong inc error poly: " << u_strong << std::endl;

    return 0;
}
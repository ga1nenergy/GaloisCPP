//
// Created by ga1nenergy on 23.07.2019.
//

#include "gfpoly.h"
#include "gfelement.h"

#include <vector>

#ifndef GALOISCPP_GFALGORITHMS_H
#define GALOISCPP_GFALGORITHMS_H

namespace galoiscpp {
    GFpoly newton_identitites(const std::vector<GFelement> &pow_sum);
    GFpoly lagrange_interpolation();
    std::vector<std::vector<Fint>> gauss_method_modulus(const std::vector<std::vector<Fint>> &matrix, Fint modulus);
}

#endif //GALOISCPP_GFALGORITHMS_H

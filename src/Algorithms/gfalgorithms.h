//
// Created by ga1nenergy on 23.07.2019.
//

#ifndef GALOISCPP_GFALGORITHMS_H
#define GALOISCPP_GFALGORITHMS_H

#include "GFpoly/gfpoly.h"
#include "GFelement/gfelement.h"
#include "typedefs.h"

#include <vector>

namespace galoiscpp {
    GFpoly newton_identitites(const std::vector<GFelement> &pow_sum);
    GFpoly lagrange_interpolation();

    std::vector<std::vector<Fint>> gauss_method_modulus(const std::vector<std::vector<Fint>> &matrix, Fint modulus);
    std::vector<std::vector<GFelement>> gauss_method(const std::vector<std::vector<GFelement>> &matrix);


}

#endif //GALOISCPP_GFALGORITHMS_H

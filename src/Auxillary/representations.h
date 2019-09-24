//
// Created by ga1nenergy on 22.07.2019.
//

#ifndef GALOISCPP_REPRESENTATIONS_H
#define GALOISCPP_REPRESENTATIONS_H

#include "GFpoly/gfpoly.h"

namespace galoiscpp {
    GFpoly roots_to_poly(const std::vector<GFelement> &rts);
    GFpoly roots_to_poly(const GaloisField *field, const std::vector<int> &rts);

    GFpoly locator_polynomial(const std::vector<GFelement> &locators, const std::vector<int> &msg);
    GFpoly locator_polynomial(const GaloisField *field, const std::vector<Fint> &locators, const std::vector<int> &msg);
}

#endif //GALOISCPP_REPRESENTATIONS_H

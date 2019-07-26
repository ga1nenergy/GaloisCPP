//
// Created by ga1nenergy on 22.07.2019.
//

#include "gfpoly.h"

#ifndef GALOISCPP_REPRESENTATIONS_H
#define GALOISCPP_REPRESENTATIONS_H

namespace galoiscpp {
    GFpoly roots_to_poly(const std::vector<GFelement> &rts);
    GFpoly roots_to_poly(GaloisField *field, const std::vector<Fint> &rts);

    GFpoly locator_polynomial(const std::vector<GFelement> &locators, const std::vector<Fint> &msg);
    GFpoly locator_polynomial(GaloisField *field, const std::vector<Fint> &locators, const std::vector<Fint> &msg);
}

#endif //GALOISCPP_REPRESENTATIONS_H

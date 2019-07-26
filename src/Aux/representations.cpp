//
// Created by ga1nenergy on 22.07.2019.
//

#include "representations.h"
#include "../../include/galoisfield.h"

namespace galoiscpp {
    GFpoly roots_to_poly(GaloisField *field, const std::vector<Fint> &rts_int) {
        std::vector<GFelement> rts;
        for (auto & elem : rts_int) rts.emplace_back(field, elem);

        GFpoly res = roots_to_poly(rts);
        return res;
    }

    /* TODO
     * 1) add roots field compatibility check
     */
    GFpoly roots_to_poly(const std::vector<GFelement> &rts) {
        GaloisField *field = rts[0].getField();

        GFpoly res(field, 0);
        res[0] = GFelement(field, 1);

        for (auto &rt : rts) {
            GFpoly root_poly(field, 1);
            if (rt == 0) {
                root_poly[1] = GFelement(field, 1);
            } else {
                root_poly[0] = GFelement(field, 1);
                root_poly[1] = -rt.inverse();
            }

            res = res * root_poly;
        }

        return res;
    }

    GFpoly locator_polynomial(const std::vector<GFelement> &locators, const std::vector<Fint> &msg) {
        if (locators.size() != msg.size())
            throw std::logic_error("Array sizes must be equal");

        std::vector<GFelement> all_rts;
        for (size_t i = 0; i < locators.size(); i++) {
            for (Fint j = 0; j < msg[i]; j++)
                all_rts.push_back(locators[i]);
        }

        GFpoly res = roots_to_poly(all_rts);
        return res;
    }

    GFpoly locator_polynomial(GaloisField *field, const std::vector<Fint> &locators, const std::vector<Fint> &msg) {
        std::vector<GFelement> locators_gf;
        for (auto & elem : locators) locators_gf.emplace_back(field, elem);

        GFpoly res = locator_polynomial(locators_gf, msg);
        return res;
    }
}
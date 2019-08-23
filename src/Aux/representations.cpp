//
// Created by ga1nenergy on 22.07.2019.
//

#include "representations.h"
#include "GaloisField/galoisfield.h"

namespace galoiscpp {
    GFpoly roots_to_poly(const GaloisField &field, const std::vector<int> &rts_int) {
        std::vector<GFelement> rts;
        for (auto & elem : rts_int) rts.emplace_back(field, elem);

        GFpoly res = roots_to_poly(rts);
        return res;
    }

    /* TODO
     * 1) add roots field compatibility check
     */
    GFpoly roots_to_poly(const std::vector<GFelement> &rts) {
        if (rts.empty()) {
            throw std::logic_error("Input array is empty");
        }

        GaloisField field = rts[0].getField();

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

    GFpoly locator_polynomial(const std::vector<GFelement> &locators, const std::vector<int> &msg) {
        if (locators.size() != msg.size())
            throw std::logic_error("Array sizes must be equal");

        if (locators.empty() || msg.empty()) {
            throw std::logic_error("Input array is empty");
        }

        std::vector<GFelement> all_rts;
        for (size_t i = 0; i < locators.size(); i++) {
            for (Fint j = 0; j < msg[i]; j++)
                all_rts.push_back(locators[i]);
        }

        GFpoly res;
        if (!all_rts.empty()) {
            res = roots_to_poly(all_rts);
        } else {
            res = GFpoly(locators[0].getField(), 0);
            res[0] = 1;
        }

        return res;
    }

    GFpoly locator_polynomial(const GaloisField &field, const std::vector<Fint> &locators, const std::vector<int> &msg) {
        std::vector<GFelement> locators_gf;
        for (auto & elem : locators) locators_gf.emplace_back(field, elem);

        GFpoly res = locator_polynomial(locators_gf, msg);
        return res;
    }
}
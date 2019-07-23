//
// Created by ga1nenergy on 19.07.2019.
//

#ifndef GALOISCPP_GFPOLY_H
#define GALOISCPP_GFPOLY_H

#include <vector>
#include "gfelement.h"

// LSB goes first!

namespace galoiscpp {
class GFpoly {
    GaloisField* field;
    std::vector<GFelement> coefs;
    Fint degree;

    GFpoly reverse() const;

public:
    GFpoly();
    GFpoly(const std::vector<GFelement>& coefs);
    GFpoly(GaloisField* gf);
    GFpoly(GaloisField* gf, const std::vector<Fint>& coefs);
    GFpoly(GaloisField *gf, Fint deg);
    GFpoly(const GFpoly &poly);
    ~GFpoly();

    GFpoly& operator=(const GFpoly& lhs);

    GFpoly operator+(const GFpoly &lhs) const;
    GFpoly operator-(const GFpoly &lhs) const;
    GFpoly operator-() const;

    friend GFpoly operator*(const GFpoly &lhs, const GFpoly &rhs);
    friend GFpoly operator*(const GFelement &lhs, const GFpoly &rhs);
    friend GFpoly operator*(const GFpoly &lhs, const GFelement &rhs);

    /* TODO
     * 1) Flip remainder (LSB goes first)
     */
    std::vector<GFpoly> operator/(const GFpoly &lhs) const;
    GFpoly operator/(const GFelement &lhs) const;

    /* TODO
     * 1) add poly field compatibility here
     */
    friend bool operator==(const GFpoly& lhs, const GFpoly& rhs);
    friend bool operator!=(const GFpoly& lhs, const GFpoly& rhs);
    friend bool operator>(const GFpoly& lhs, const GFpoly& rhs);
    friend bool operator<(const GFpoly& lhs, const GFpoly& rhs);
    friend bool operator>=(const GFpoly& lhs, const GFpoly& rhs);
    friend bool operator<=(const GFpoly& lhs, const GFpoly& rhs);

    friend std::ostream &operator<<(std::ostream &os, const GFpoly &poly);

    GFelement& operator[](int idx);
    GFelement operator[](int idx) const;

    std::vector<GFelement> roots() const;
    static std::vector<GFelement> roots(const GFpoly& poly);

    GFelement polyval(const GFelement& elem) const;
    GFelement polyval(Fint elem) const;
    static GFelement polyval(const GFpoly& poly, const GFelement& elem);
    static GFelement polyval(const GFpoly& poly, Fint elem);

    /* TODO
     * 1) check euclid methods
     */
    static GFpoly euclid(const GFpoly &op1, const GFpoly &op2, Fint degree);
    static std::pair<GFelement, GFpoly> extended_euclid(const GFpoly &op1, const GFpoly &op2);

    GFpoly inverse(const GFpoly &mod) const;
    static GFpoly inverse(const GFpoly &poly, const GFpoly &mod);

    bool empty();


    /* TODO:
     * 1) inverse           DONE
     * 2) roots             DONE
     * 3) polyval           DONE
     * 4) euclid            DONE
     * 5) extended_euclid   DONE
     */

    // Getters
    GaloisField* getField() const;
};
}

#endif //GALOISCPP_GFPOLY_H

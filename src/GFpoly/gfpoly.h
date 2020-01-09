//
// Created by ga1nenergy on 19.07.2019.
//

#ifndef GALOISCPP_GFPOLY_H
#define GALOISCPP_GFPOLY_H

#include <vector>
#include "GaloisField/galois_field_prime.h"
#include "GFelement/gfelement.h"
#include "Auxillary/structures.h"

// LSB goes first!

namespace galoiscpp {
class GFpoly {
    const GaloisField* field;
    std::vector<GFelement> coefs;
    Fint degree;

    GFpoly reverse() const;
    GFpoly cut_zeros() const;

public:
    GFpoly();
    GFpoly(const std::vector<GFelement>& coefs);
    GFpoly(const GaloisField *gf);
    GFpoly(const GaloisField *gf, const std::vector<Fint>& coefs);
    GFpoly(const GaloisField *gf, Fint deg);
    GFpoly(const GFpoly &poly);
    ~GFpoly();

    GFpoly& operator=(const GFpoly& rhs);

    GFpoly operator+(const GFpoly &rhs) const;
    GFpoly operator+(const GFelement &rhs) const;
    GFpoly operator+(Fint rhs) const;
    GFpoly operator+() const;
    GFpoly operator-(const GFpoly &rhs) const;
    GFpoly operator-(const GFelement &rhs) const;
    GFpoly operator-(Fint rhs) const;
    GFpoly operator-() const;

    friend GFpoly operator*(const GFpoly &lhs, const GFpoly &rhs);
    friend GFpoly operator*(const GFelement &lhs, const GFpoly &rhs);
    friend GFpoly operator*(const GFpoly &lhs, const GFelement &rhs);

    /* TODO
     * 1) Flip remainder (LSB goes first)
     */
    std::vector<GFpoly> operator/(const GFpoly &rhs) const;
    GFpoly operator/(const GFelement &rhs) const;

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

    static GFpoly conv(const GFpoly &op1, const GFpoly &op2);
    static GFpoly convmod(const GFpoly &op1, const GFpoly &op2, const GFpoly &mod);
    static std::vector<GFpoly> deconv(const GFpoly &op1, const GFpoly &op2);
    static GFpoly deconvq(const GFpoly &op1, const GFpoly &op2);
    static GFpoly deconvr(const GFpoly &op1, const GFpoly &op2);
    static GFpoly deconvmod(const GFpoly &op1, const GFpoly &op2, const GFpoly &mod);

    GFpoly deriv() const;
    static GFpoly deriv (const GFpoly &poly);

    static GFpoly trace_as_poly(const GaloisField *field);

    std::vector<GFpoly> find_splitting_functions() const;
    static std::vector<GFpoly> find_splitting_functions(const GFpoly &poly);

    std::vector<GFelement> roots() const;
    static std::vector<GFelement> roots(const GFpoly& poly);
    static std::vector<GFelement> roots(const std::vector<Fint> &coefs, const std::vector<GFelement> &set);

    std::vector<GFelement> roots_exhaustive() const;
    static std::vector<GFelement> roots_exhaustive(const GFpoly &poly);
    std::vector<GFelement> roots_exhaustive(const std::vector<GFelement> &set) const;
    static std::vector<GFelement> roots_exhaustive(const GFpoly &poly, const std::vector<GFelement> &set);
    std::vector<GFelement> roots_prime() const;
    static std::vector<GFelement> roots_prime(const GFpoly &poly);
    std::vector<GFelement> roots1() const;
    static std::vector<GFelement> roots1(const GFpoly& poly);
    std::vector<GFelement> roots2() const;
    static std::vector<GFelement> roots2(const GFpoly& poly);

    GFelement polyval(const GFelement& elem) const;
    GFelement polyval(Fint elem) const;
    static GFelement polyval(const GFpoly& poly, const GFelement& elem);
    static GFelement polyval(const GFpoly& poly, Fint elem);
    static GFelement polyval(const std::vector<Fint> &coefs, const GFelement &elem);

    static GFpoly euclid(const GFpoly &op1, const GFpoly &op2, Fint degree);
    static std::vector<GFpoly> euclid_list(const GFpoly &op1, const GFpoly &op2);
    static std::pair<GFelement, GFpoly> extended_euclid(const GFpoly &op1, const GFpoly &op2);
    static GFpoly gcd(const GFpoly &op1, const GFpoly &op2);

    static std::vector<std::vector<GFelement>> get_cosets(GaloisField *field, std::vector<int> *elem_idx = nullptr);
    static std::vector<GFpoly> get_primitive_polys(GaloisField *field, GaloisFieldPrime *prime_field);

    GFpoly inverse(const GFpoly &mod) const;
    static GFpoly inverse(const GFpoly &poly, const GFpoly &mod);

    bool empty();
    bool isZero();


    /* TODO:
     * 1) inverse           DONE
     * 2) roots             DONE
     * 3) polyval           DONE
     * 4) euclid            DONE
     * 5) extended_euclid   DONE
     * 6) roots_alg
     * 7) gcd
     * 8) trace_as_poly     DONE
     */

    // Getters
    const GaloisField* getField() const;
    Fint getDegree() const;
    std::vector<GFelement> getCoefs() const;
    std::vector<Fint> getCoefsRaw() const;
};
}

#endif //GALOISCPP_GFPOLY_H

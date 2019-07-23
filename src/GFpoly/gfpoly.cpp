//
// Created by ga1nenergy on 19.07.2019.
//

#include "gfpoly.h"
#include <vector>
#include <algorithm>

namespace galoiscpp {
    GFpoly::GFpoly() = default;

    GFpoly::GFpoly(const vector<GFelement> &coefs) {
        this->field = coefs[0].getField();
        for (auto & coef : coefs) {
            if (coef.getField() != field)
                throw std::logic_error("Fields mismatch");
            this->coefs.push_back(coef);
        }
        this->degree = static_cast<Fint>(coefs.size()) - 1;
    }

    GFpoly::GFpoly(GaloisField* field) {
        this->field= field;
    }

    GFpoly::GFpoly(GaloisField* field, const vector<Fint> &coefs) {
        this->field = field;
        this->degree = static_cast<Fint>(coefs.size()) - 1;
        for (auto & coef : coefs) this->coefs.push_back(GFelement(field, coef));
    }

    GFpoly::GFpoly(GaloisField *field, Fint degree) {
        this->field = field;
        this->degree = degree;
        for (int i = 0; i < degree + 1; i++) coefs.push_back(GFelement(field, 0));
    }

    GFpoly::GFpoly(const GFpoly &poly) = default;

    GFpoly::~GFpoly() = default;

    GFpoly GFpoly::reverse() const {
        std::vector<GFelement> reversed_coef = this->coefs; std::reverse(reversed_coef.begin(), reversed_coef.end());
        GFpoly res(reversed_coef);
        return res;
    }


    GFpoly& GFpoly::operator=(const GFpoly &lhs) {
        field = lhs.field;
        coefs = lhs.coefs;
        degree = lhs.degree;

        return *this;
    }


    GFpoly GFpoly::operator+(const GFpoly &lhs) const {
        GFpoly res, op;
        (this->degree > lhs.degree) ? ({res = *this; op = lhs;}) : ({res = lhs; op = *this;});

        for (int i = 0; i <= op.degree; i++) {
            res[i] = res[i] + op[i];
        }

        while (res.coefs[res.degree] == 0 && res.degree != 0) {
            res.coefs.pop_back();
            res.degree--;
        }

        return res;
    }

    GFpoly GFpoly::operator-(const GFpoly &lhs) const {
        GFpoly res, op;
        (this->degree > lhs.degree) ? ({res = *this; op = lhs;}) : ({res = lhs; op = *this;});

        for (int i = 0; i <= op.degree; i++) {
            res[i] = res[i] - op[i];
        }

        while (res.coefs[res.degree] == 0 && res.degree != 0) {
            res.coefs.pop_back();
            res.degree--;
        }

        return res;
    }

    GFpoly GFpoly::operator-() const {
        GFpoly res = *this;

        for (auto &coef : res.coefs) {
            coef = -coef;
        }

        return res;
    }

    GFpoly operator*(const GFpoly &lhs, const GFpoly &rhs) {
        GFpoly res(rhs.getField(), rhs.degree + lhs.degree);

        GFpoly op1 = rhs, op2 = lhs;

        for (size_t i = 0; i < op1.coefs.size(); i++) {
            for (size_t j = 0; j < op2.coefs.size(); j++) {
                res[i+j] = res[i+j] + op1[i] * op2[j];
            }
        }

        return res;
    }

    GFpoly operator*(const GFelement &lhs, const GFpoly &rhs) {
        GFpoly res = rhs;

        for (auto & elem : res.coefs)
            elem = elem * lhs;
        return res;
    }

    GFpoly operator*(const GFpoly &lhs, const GFelement &rhs) {
        return rhs * lhs;
    }

    /* TODO:
     * 1) Add poly field compatibility check
     * 2)
     */
    std::vector<GFpoly> GFpoly::operator/(const GFpoly &rhs) const {
        std::vector<GFpoly> qr(2);
        GFpoly r = *this;

        GFpoly denom = rhs;
        std::vector<GFelement> q;

        if (this->degree < rhs.degree) {
            qr[0] = GFpoly(field, 0);
            qr[1] = *this;
            return qr;
        }

        while (r.degree >= denom.degree) {
            GFelement coef = r[r.degree] * denom[denom.degree].inverse();
            GFpoly op = coef * denom;

            for (int i = 0; i < op.degree + 1; i++) {
                r[i + r.degree - op.degree] = r[i + r.degree - op.degree] - op[i];
            }
            r.coefs.pop_back();
            r.degree--;

            q.insert(q.begin(), coef);
        }

        while (r[r.degree] == 0 && r.degree != 0) {
            r.coefs.pop_back();
            r.degree--;
        }

        qr[0] = GFpoly(q);
        qr[1] = r;

        return qr;
    }

    GFpoly GFpoly::operator/(const GFelement &lhs) const {
        GFpoly res = *this;

        for (auto & elem : res.coefs)
            elem = elem / lhs;
        return res;
    }

    GFelement& GFpoly::operator[](int idx) {
        return coefs[idx];
    }

    GFelement GFpoly::operator[](int idx) const {
        return coefs[idx];
    }

    // Getters
    GaloisField* GFpoly::getField() const {
        if (!coefs.empty()) {
            return coefs[0].getField();
        } else {
            throw std::logic_error("Poly is empty");
        }
    }

    bool operator==(const GFpoly &lhs, const GFpoly &rhs) {
        if (lhs.degree != rhs.degree)
            return false;

        for (int i = 0; i < lhs.degree + 1; i++)
            if (lhs[i] != rhs[i])
                return false;

        return true;
    }

    bool operator!=(const GFpoly &lhs, const GFpoly &rhs) {
        return !(lhs == rhs);
    }

    bool operator>(const GFpoly &lhs, const GFpoly &rhs) {
        if (lhs.degree > rhs.degree)
            return true;

        for (int i = 0; i < lhs.degree + 1; i++)
            if (lhs[i] <= rhs[i])
                return false;

        return true;
    }

    bool operator<(const GFpoly &lhs, const GFpoly &rhs) {
        if (lhs.degree < rhs.degree)
            return true;

        for (int i = 0; i < lhs.degree + 1; i++)
            if (lhs[i] >= rhs[i])
                return false;

        return true;
    }

    bool operator>=(const GFpoly &lhs, const GFpoly &rhs) {
        return !(lhs < rhs);
    }

    bool operator<=(const GFpoly &lhs, const GFpoly &rhs) {
        return !(lhs > rhs);
    }

    std::ostream &operator<<(std::ostream &os, const GFpoly &poly) {
        for (auto const & elem : poly.coefs)
            os << elem << " ";
        return os;
    }

    std::vector<GFelement> GFpoly::roots() const {
        std::vector<GFelement> rts;

        for (int i = 0; i < field->get_size(); i++) {
            if (this->polyval(i) == 0) {
                GFelement rt = GFelement(field, i);

                GFpoly q = *this, r;
                while (true) {
                    GFpoly root_poly(field, 1);
                    if (i == 0) {
                        root_poly[0] = GFelement(field, 0);
                        root_poly[1] = GFelement(field, 1);
                    } else {
                        root_poly[0] = GFelement(field, 1);
                        root_poly[1] = -rt.inverse();
                    }

                    auto qr = q / root_poly;
                    q = qr[0]; cout << "q: " << q << endl;
                    r = qr[1]; cout << "r: " << r << endl;

                    if (r[0] != 0 || r.degree != 0)
                        break;

                    rts.push_back(rt);
                }
            }
        }

        return rts;
    }

    std::vector<GFelement> GFpoly::roots(const GFpoly &poly) {
        return poly.roots();
    }

    GFelement GFpoly::polyval(Fint elem) const {
        return this->polyval(GFelement(field, elem));
    }

    GFelement GFpoly::polyval(const GFelement& elem) const {
        GFelement s(coefs[0]);

        for (int i = 1; i < degree + 1; i++) {
            s = s + coefs[i] * elem.power(i);
        }
        return s;
    }

    bool GFpoly::empty() {
        return coefs.empty();
    }

    GFelement GFpoly::polyval(const GFpoly &poly, const GFelement &elem) {
        return poly.polyval(elem);
    }

    GFelement GFpoly::polyval(const GFpoly &poly, Fint elem) {
        return poly.polyval(elem);
    }

    GFpoly GFpoly::euclid(const GFpoly &op1, const GFpoly &op2, Fint degree) {
        GFpoly num, denom;
        (op1.degree > op2.degree) ? ({num = op1; denom = op2;}) : ({num = op2; denom = op1;});

        while (denom.degree >= degree) {
            auto qr = num / denom;
            num = denom;
            denom = qr[1];
        }

        return denom;
    }

    /* TODO
     * 1) add field compatibility check
     */

    std::pair<GFelement, GFpoly> GFpoly::extended_euclid(const GFpoly &op1, const GFpoly &op2) {
        GaloisField *field = op1.getField();
        GFpoly r2 = op2, r1 = op1;
        GFpoly y2(field, 0); y2[0] = GFelement(field, 0);
        GFpoly y1(field, 0); y1[0] = GFelement(field, 1);

        int i = 1;
        cout << "Step " << i << endl;
        cout << "r2: " << r2 << endl;
        cout << "r1: " << r1 << endl;
        cout << "y2: " << y2 << endl;
        cout << "y1: " << y1 << endl;
        cout << endl;

        while (r1.degree != 0) {
            std::vector<GFpoly> qr = r2 / r1;

            GFpoly buf = y1;
            y1 = y2 - qr[0] * y1;
            y2 = buf;
            r2 = r1;
            r1 = qr[1];

            i++;
            cout << "Step " << i << endl;
            cout << "r2: " << r2 << endl;
            cout << "r1: " << r1 << endl;
            cout << "y2: " << y2 << endl;
            cout << "y1: " << y1 << endl;
        }

        std::pair<GFelement, GFpoly> res = {r1[0], y1};
        return res;
    }

    GFpoly GFpoly::inverse(const GFpoly &mod) const {
        std::pair<GFelement, GFpoly> p = extended_euclid(*this, mod);
        GFpoly res = std::get<0>(p).inverse() * std::get<1>(p);

        return res;
    }

    GFpoly GFpoly::inverse(const GFpoly &poly, const GFpoly &mod) {
        GFpoly res = poly.inverse(mod);
        return res;
    }

}

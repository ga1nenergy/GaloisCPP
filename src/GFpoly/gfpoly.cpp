//
// Created by ga1nenergy on 19.07.2019.
//

#include "gfpoly.h"
#include "Auxillary/structures.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>
#include <galoiscpp.h>

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

    GFpoly::GFpoly(const GaloisField *field) {
        this->field = field;
    }

    GFpoly::GFpoly(const GaloisField *field, const vector<Fint> &coefs) {
        this->field = field;
        this->degree = static_cast<Fint>(coefs.size()) - 1;
        for (auto & coef : coefs) this->coefs.emplace_back(field, coef);
    }

    GFpoly::GFpoly(const GaloisField *field, Fint degree) {
        this->field = field;
        this->degree = degree;
        for (int i = 0; i < degree + 1; i++) coefs.emplace_back(field, 0);
    }

    GFpoly::GFpoly(const GFpoly &poly) = default;

    GFpoly::~GFpoly() = default;

    GFpoly GFpoly::reverse() const {
        std::vector<GFelement> reversed_coef = this->coefs; std::reverse(reversed_coef.begin(), reversed_coef.end());
        GFpoly res(reversed_coef);
        return res;
    }

    GFpoly GFpoly::cut_zeros() const {
        auto i = degree;
        while (coefs[i] == 0) i--;

        auto res = *this;
        res.coefs.resize(i + 1);
        return res;
    }

    GFpoly& GFpoly::operator=(const GFpoly &rhs) {
        field = rhs.field;
        coefs = rhs.coefs;
        degree = rhs.degree;

        return *this;
    }


    GFpoly GFpoly::operator+(const GFpoly &rhs) const {
        GFpoly res, op;
        (this->degree > rhs.degree) ? ({res = *this; op = rhs;}) : ({res = rhs; op = *this;});

        for (int i = 0; i <= op.degree; i++) {
            res[i] = res[i] + op[i];
        }

        while (res.coefs[res.degree] == 0 && res.degree != 0) {
            res.coefs.pop_back();
            res.degree--;
        }

        return res;
    }

    GFpoly GFpoly::operator+(const GFelement &rhs) const {
        auto res = *this;

        if (this->getDegree() >= 0) {
            res[0] = res[0] + rhs;
        } else {
            throw std::logic_error("Poly is uninitialized");
        }

        return res;
    }

    GFpoly GFpoly::operator+(Fint rhs) const {
        auto res = *this;

        if (this->getDegree() >= 0) {
            res[0] = res[0] + rhs;
        } else {
            throw std::logic_error("Poly is uninitialized");
        }

        return res;
    }

    GFpoly GFpoly::operator+() const {
        return *this;
    }

    GFpoly GFpoly::operator-(const GFpoly &rhs) const {
        GFpoly res;

        if (this->degree > rhs.degree) {
            res = *this;
            for (int i = 0; i <= rhs.degree; i++) {
                res[i] = res[i] - rhs[i];
            }
        } else {
            res = rhs;
            for (int i = 0; i <= this->degree; i++) {
                res[i] = this->coefs[i] - res[i];
            }
            for (int i = this->degree + 1; i <= res.degree; i++)
                res[i] = -res[i];
        }

        while (res.coefs[res.degree] == 0 && res.degree != 0) {
            res.coefs.pop_back();
            res.degree--;
        }

        return res;
    }

    GFpoly GFpoly::operator-(const GFelement &rhs) const {
        auto res = *this;

        if (this->getDegree() >= 0) {
            res[0] = res[0] - rhs;
        } else {
            throw std::logic_error("Poly is uninitialized");
        }

        return res;
    }

    GFpoly GFpoly::operator-(Fint rhs) const {
        auto res = *this;

        if (this->getDegree() >= 0) {
            res[0] = res[0] - rhs;
        } else {
            throw std::logic_error("Poly is uninitialized");
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

//        cout << "op1: " << op1 << endl;
//        cout << "op2: " << op2 << endl;

        for (size_t i = 0; i < op1.coefs.size(); i++) {
            for (size_t j = 0; j < op2.coefs.size(); j++) {
                res[i+j] = res[i+j] + op1[i] * op2[j];
//                std::cout << res << ", " << op1[i] * op2[j] << std::endl;
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
//            std::cout << "coef: " << coef << std::endl;
            GFpoly op = coef * denom;
//            std::cout << "op: " << op << std::endl;
//            std::cout << "r before: " << r << std::endl;

//            GFelement op1(field, 0);
//            GFelement op2(field, 6);
//
//            op1 - op2;

            for (int i = op.degree; i >= 0; i--) {
//                cout << r[i + r.degree - op.degree] << " " << op[i] << endl;
                r[i + r.degree - op.degree] = r[i + r.degree - op.degree] - op[i];
//                cout << r[i + r.degree - op.degree] << endl;
            }
            r.coefs.pop_back();
            r.degree--;
//            std::cout << "r after: " << r << std::endl;

            q.insert(q.begin(), coef);
//            std::cout << "q_iter: ";
//            std::copy(q.begin(), q.end(), std::ostream_iterator<GFelement>(std::cout, " "));
//            std::cout << std::endl;
        }

        while (r[r.degree] == 0 && r.degree != 0) {
            r.coefs.pop_back();
            r.degree--;
        }

//        auto i = std::find_if(std::rbegin(r.coefs), std::rend(r.coefs), [](auto& v) { return v != 0; } );
//        auto diff = std::distance(std::begin(r.coefs), i.base());
//
//        if (diff == 0){
//            r.coefs.resize(diff + 1);
//            r.degree = diff;
//        } else {
//            r.coefs.resize(diff);
//            r.degree = diff - 1;
//        }
//
//        i = std::find_if(std::rbegin(q), std::rend(q), [](auto& v) { return v != 0; } );
//        diff = std::distance(std::begin(q), i.base());
//        (diff == 0) ? (q.resize(diff + 1)) : (q.resize(diff));

        qr[0] = GFpoly(q);
        qr[1] = r;

        return qr;
    }

    GFpoly GFpoly::operator/(const GFelement &rhs) const {
        GFpoly res = *this;

        for (auto & elem : res.coefs)
            elem = elem / rhs;
        return res;
    }

    GFelement& GFpoly::operator[](int idx) {
        return coefs[idx];
    }

    GFelement GFpoly::operator[](int idx) const {
        return coefs[idx];
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

    GFpoly GFpoly::conv(const GFpoly &op1, const GFpoly &op2) {
        return op1 * op2;
    }

    GFpoly GFpoly::convmod(const GFpoly &op1, const GFpoly &op2, const GFpoly &mod) {
        auto qr = (op1 * op2) / mod;
        return qr[1];
    }

    std::vector<GFpoly> GFpoly::deconv(const GFpoly &op1, const GFpoly &op2) {
        return op1 / op2;
    }

    GFpoly GFpoly::deconvq(const GFpoly &op1, const GFpoly &op2) {
        auto qr = op1 / op2;
        return qr[0];
    }

    GFpoly GFpoly::deconvr(const GFpoly &op1, const GFpoly &op2) {
        auto qr = op1 / op2;
        return qr[1];
    }

    GFpoly GFpoly::deconvmod(const GFpoly &op1, const GFpoly &op2, const GFpoly &mod) {
        auto qr = (op1 * op2.inverse(mod)) / mod;
        return qr[1];
    }

    GFpoly GFpoly::deriv() const {
        if (degree == 0) {
            GFpoly res(field, 0); res[0] = 0;
            return res;
        }

        GFpoly res(field, degree - 1);
        for (auto i = 0; i <= res.getDegree(); i++) {
            res[i] = (*this)[i+1].sum_times(i+1);
        }

        return res.cut_zeros();
    }

    GFpoly GFpoly::deriv(const GFpoly &poly) {
        return poly.deriv();
    }

    GFpoly GFpoly::trace_as_poly(const GaloisField *field) {
        GFpoly tr(field, pow(field->getModulus(), field->getDimension() - 1));

        for (auto i = 0; i < field->getDimension(); i++) {
            tr[(int)pow(field->getModulus(), i)] = 1;
        }

        return tr;
    }

    // a shortened implementation of Berlekamp algorithm
    // TODO: finish function
    std::vector<GFpoly> GFpoly::find_splitting_functions() const {
        // Find all x^(iq) mod f(x), i=0,1,...,n-1 and put them as rows of matrix B
        std::vector<std::vector<GFelement>> B(degree);
        std::vector<GFpoly> factors(degree);
        for (auto i = 0; i < degree; i++) {
            GFpoly px(field, i * field->get_size()); px[px.degree] = 1;
            px = deconvr(px, *this);
            factors[i] = px;

            std::vector<GFelement> row(degree);
            std::copy(px.coefs.begin(), px.coefs.end(), row.begin());
            std::for_each(row.begin() + px.coefs.size(), row.end(), [&](auto &e){e = GFelement(field, 0);});
            B[i] = row;
        }

        // Find rank of B and splitting functions
        /*auto solution = gauss_method(B);
        std::vector<GFelement> zero_vector(degree);
        std::for_each(zero_vector.begin(), zero_vector.end(), [&](auto &e){e = GFelement(field, 0);});

        for (auto i = solution.begin(); i < solution.end(); i++) {
            if (*i == zero_vector)
        }
        auto rank = std::min(solution.size(), solution[0].size());*/
        {};
    }

    std::vector<GFpoly> GFpoly::find_splitting_functions(const GFpoly &poly) {
        return std::vector<GFpoly>();
    }

    // TODO: implement adaptive choise of roots finding algorithm
    std::vector<GFelement> GFpoly::roots() const {
        return roots_exhaustive();
    }

    std::vector<GFelement> GFpoly::roots(const GFpoly &poly) {
        return poly.roots();
    }

    std::vector<GFelement> GFpoly::roots(const std::vector<Fint> &coefs, const std::vector<GFelement> &set) {
        std::vector<GFelement> rts;

        for (auto const &rt : set) {
            if (polyval(coefs, rt) == 0) {
                rts.push_back(rt);
            }
        }

        return rts;
    }

    std::vector<GFelement> GFpoly::roots_exhaustive() const {
        std::vector<GFelement> rts;

        for (int i = 0; i < field->get_size(); i++) {
            if (this->polyval(i) == 0) {
                GFelement rt = GFelement(field, i);

                GFpoly root_poly(field, 1);
                if (i == 0) {
                    root_poly[0] = GFelement(field, 0);
                    root_poly[1] = GFelement(field, 1);
                } else {
                    root_poly[0] = GFelement(field, 1);
                    root_poly[1] = -rt.inverse();
                }

//                std::cout << "roots poly: " << root_poly << std::endl;

                GFpoly q = *this, r;
                while (true) {
//                    std::cout << "q: " << q << std::endl;
                    auto qr = q / root_poly;
                    q = qr[0]; //cout << "q: " << q << endl;
                    r = qr[1]; //cout << "r: " << r << endl;

                    if (r[0] != 0 || r.degree != 0)
                        break;

                    rts.push_back(rt);
                }
            }
        }

        return rts;
    }

    std::vector<GFelement> GFpoly::roots_exhaustive(const std::vector<GFelement> &set) const {
        std::vector<GFelement> rts;

        for (auto const &rt : set) {
            if (this->polyval(rt) == 0) {
                GFpoly root_poly(field, 1);
                if (rt == 0) {
                    root_poly[0] = GFelement(field, 0);
                    root_poly[1] = GFelement(field, 1);
                } else {
                    root_poly[0] = GFelement(field, 1);
                    root_poly[1] = -rt.inverse();
                }

//                std::cout << "roots poly: " << root_poly << std::endl;

                GFpoly q = *this, r;
                while (true) {
                    auto qr = q / root_poly;
                    q = qr[0]; // cout << "q: " << q << endl;
                    r = qr[1]; // cout << "r: " << r << endl;

                    if (r[0] != 0 || r.degree != 0 || q[q.getDegree()] == 0)
                        break;

                    rts.push_back(rt);
                }
            }
        }

        return rts;
    }

    std::vector<GFelement> GFpoly::roots_exhaustive(const GFpoly &poly, const std::vector<GFelement> &set) {
        return poly.roots_exhaustive(set);
    }

    std::vector<GFelement> GFpoly::roots_exhaustive(const GFpoly &poly) {
        return poly.roots_exhaustive();
    }

    std::vector<GFelement> GFpoly::roots_prime() const {
        return std::vector<GFelement>();
    }

    std::vector<GFelement> GFpoly::roots_prime(const GFpoly &poly) {
        return poly.roots_prime();
    }

    // suitable for cases when p is small and m is large
    // nop. exhaustive search: 2^m, algorithm: 2^(deg - 1)
    // TODO; finish function
    std::vector<GFelement> GFpoly::roots1() const {
        auto tr_poly = trace_as_poly(this->field);

        // main_poly = x^q + x
        GFpoly main_poly(field, field->get_size());
        main_poly[1] = 1; main_poly[main_poly.getDegree()] = 1;

        {};
    }

    std::vector<GFelement> GFpoly::roots1(const GFpoly &poly) {
        return poly.roots1();
    }

    // suitable for cases when p^m is large.
    // TODO: finish function
    std::vector<GFelement> GFpoly::roots2() const {
        std::vector<GFpoly> f_k_array;
        GFpoly F(field, 0); F[0] = 1;
        for (auto k = 0; k < field->getDimension(); k++) {
            auto f_k = *this;
            for (auto j = 0; j <= f_k.getDegree(); j++) {
                f_k[j] = f_k[j].power(pow(field->getModulus(), k));
            }
            F = F * f_k;
        }
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

    GFelement GFpoly::polyval(const std::vector<Fint> &coefs, const GFelement &elem) {
        auto field = elem.getField();
        auto degree = coefs.size() - 1;

        GFelement s(field, coefs[0]);

        for (int i = 1; i < degree + 1; i++) {
            s = s + elem.power(i).sum_times(coefs[i]);
        }
        return s;
    }

    bool GFpoly::empty() {
        return coefs.empty();
    }

    bool GFpoly::isZero() {
        return (degree == 0 && coefs[0] == 0);
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

//        std::cout << "degree " << denom.degree << ". " << denom << std::endl;

        while (denom.degree > degree) {
            auto qr = num / denom;
            num = denom;
            denom = qr[1];

//            std::cout << "degree " << denom.degree << ". " << denom << std::endl;
        }

        return denom;
    }

    std::vector<GFpoly> GFpoly::euclid_list(const GFpoly &op1, const GFpoly &op2) {
        GFpoly num, denom;
        (op1.degree > op2.degree) ? ({num = op1; denom = op2;}) : ({num = op2; denom = op1;});

//        std::cout << "degree " << denom.degree << ". " << denom << std::endl;
        std::vector<GFpoly> rems;

        while (denom.degree > 0) {
            rems.push_back(denom);
            auto qr = num / denom;
            num = denom;
            denom = qr[1];
//            std::cout << "degree " << denom.degree << ". " << denom << std::endl;
        }

        rems.push_back(denom);

        return rems;
    }

    /* TODO
     * 1) add field compatibility check
     */

    std::pair<GFelement, GFpoly> GFpoly::extended_euclid(const GFpoly &op1, const GFpoly &op2) {
        const GaloisField *field = op1.getField();
        GFpoly r2 = op2, r1 = op1;
        GFpoly y2(field, 0); y2[0] = GFelement(field, 0);
        GFpoly y1(field, 0); y1[0] = GFelement(field, 1);

        int i = 1;
//        cout << "Step " << i << endl;
//        cout << "r2: " << r2 << endl;
//        cout << "r1: " << r1 << endl;
//        cout << "y2: " << y2 << endl;
//        cout << "y1: " << y1 << endl;


        while (r1.degree != 0) {
            std::vector<GFpoly> qr = r2 / r1;

            GFpoly buf = y1;
            y1 = y2 - qr[0] * y1;
            y2 = buf;
            r2 = r1;
            r1 = qr[1];

            i++;

//            cout << "Step " << i << endl;
//            cout << "r2: " << r2 << endl;
//            cout << "r1: " << r1 << endl;
//            cout << "y2: " << y2 << endl;
//            cout << "y1: " << y1 << endl;
        }

        std::pair<GFelement, GFpoly> res = {r1[0], y1};
        return res;
    }

    // TODO
    GFpoly GFpoly::gcd(const GFpoly &op1, const GFpoly &op2) {
        GFpoly num, denom;
        (op1.degree > op2.degree) ? ({num = op1; denom = op2;}) : ({num = op2; denom = op1;});

//        std::cout << "degree " << denom.degree << ". " << denom << std::endl;

        while (denom[0] != 0) {
            auto qr = num / denom;
            num = denom;
            denom = qr[1];

//            std::cout << "degree " << denom.degree << ". " << denom << std::endl;
        }

        return num;
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

    // TODO: make elem_idx more С++ way
    std::vector<std::vector<GFelement>> GFpoly::get_cosets(GaloisField *field, std::vector<int> *elem_idx) {
        bool LOCAL_DEBUG = false;

        std::vector<std::vector<GFelement>> cosets;
        std::vector<int> free_locators(field->get_size());
        std::iota(free_locators.begin(), free_locators.end(), 0);

        auto coset_num = 0;
        for (auto i = free_locators.begin(); i != free_locators.end(); i = free_locators.begin()) {
            std::vector<GFelement> coset;
            GFelement elem(field, *i);

            auto p_int = field->getModulus();
            auto p_power_int = p_int;
            auto current_elem = elem;
            auto j = 1;
            while (std::find(coset.begin(), coset.end(), current_elem) == coset.end()) {
                coset.push_back(current_elem);
                free_locators.erase(std::find(free_locators.begin(), free_locators.end(), current_elem.getDegree()));

                if (LOCAL_DEBUG) {
                    std::cout << "Free locators: ";
                    std::copy(free_locators.begin(), free_locators.end(), std::ostream_iterator<int>(std::cout, " "));
                    std::cout << std::endl;
                }

                if (elem_idx != nullptr) {
                    elem_idx->operator[](current_elem.getDegree()) = coset_num;
                }

                current_elem = elem.power(pow(p_int, j));
                if (LOCAL_DEBUG) {
                    std::cout << "Current power: " << pow(p_int, j) << std::endl;
                    std::cout << "Current elem: " << current_elem << std::endl;
                }

                p_power_int *= p_int;
                j++;
            }

            if (LOCAL_DEBUG) {
                std::cout << "Current coset elements: ";
                std::copy(coset.begin(), coset.end(), ostream_iterator<GFelement>(std::cout, " "));
                std::cout << std::endl;
            }

            cosets.push_back(coset);
            coset_num++;
        }

        return cosets;
    }

    std::vector<GFpoly> GFpoly::get_primitive_polys(GaloisField *field, GaloisFieldPrime *prime_field) {
        bool LOCAL_DEBUG = false;
        auto cosets = get_cosets(field);

        std::vector<GFpoly> primitive_polys;
        for (auto const &c : cosets) {
            GFpoly primitive_poly(field, 0);
            primitive_poly[0] = 1;

            if (LOCAL_DEBUG) {
                std::cout << "Current coset: ";
                std::copy(c.begin(), c.end(), std::ostream_iterator<GFelement>(std::cout, " "));
                std::cout << std::endl;
            }

            for (auto const &e : c) {
                GFpoly mul(field, 1);
                mul[0] = -e; mul[1] = 1;
                primitive_poly = primitive_poly * mul;
            }

            if (LOCAL_DEBUG) {
                std::cout << "Current primitive poly: " << primitive_poly << std::endl;
                auto rts = primitive_poly.roots();
                std::cout << "Roots: ";
                std::copy(rts.begin(), rts.end(), std::ostream_iterator<GFelement>(std::cout, " "));
                std::cout << std::endl;
            }

            for (auto i = 0; i <= primitive_poly.getDegree(); i++) { // TODO: !!!
                primitive_poly[i] = GFelement(prime_field,
                        (field->field_table_ptr()->operator[](primitive_poly[i].getDegree())[0]));
            }

            if (LOCAL_DEBUG) {
                std::cout << "New primitive poly: " << primitive_poly << std::endl;
            }

            primitive_polys.push_back(primitive_poly);
        }

        return primitive_polys;
    }

    // Getters
    const GaloisField* GFpoly::getField() const {
        if (!coefs.empty()) {
            return coefs[0].getField();
        } else {
            throw std::logic_error("Poly is empty");
        }
    }

    Fint GFpoly::getDegree() const {
        return degree;
    }

    std::vector<GFelement> GFpoly::getCoefs() const {
        return coefs;
    }

    std::vector<Fint> GFpoly::getCoefsRaw() const {
        std::vector<Fint> raw(coefs.size());
        for (auto i = 0; i < coefs.size(); i++) {
            raw[i] = coefs[i].getDegree();
        }

        return raw;
    }
}

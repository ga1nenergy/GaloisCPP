/**********************************************************************
   Project: C++ Library for General Galois Field Arithmetic

   Language: C++ 2007	   
   Author: Saied H. Khayat
   Date:   Feb 2013
   URL: https://github.com/saiedhk
   
   Copyright Notice: Free use of this library is permitted under the
   guidelines and in accordance with the MIT License (MIT).
   http://opensource.org/licenses/MIT

**********************************************************************/

#include <iostream>
#include <cmath>
#include <cstring>

using namespace std;

#include "gfelement.h"
#include "GaloisField/galois_field_prime.h"

namespace galoiscpp
{

//-----------------------------------------------------------------------------
// Operators
//-----------------------------------------------------------------------------

GFelement& GFelement::operator=(const GFelement& rhs)
{
    if (this != &rhs)
    {
        field = rhs.field;
        degree = rhs.degree;
    }
    return *this;
}


GFelement& GFelement::operator=(Fint rhs)
{
    if (field->getModulus() != -1) {
        if (rhs < field->get_size()) {
            degree = rhs;
        } else {
            throw std::logic_error("Element degree must not exceed p^m-1");
        }
    } else {
        throw std::logic_error("Assigning elem degree before initializing field");
    }
    return *this;
}



//------------------------------------------------------
GFelement GFelement::operator-() const
{
    GFelement result(field, field->subtract(0, degree));
    return result;
}



//------------------------------------------------------
bool operator==(const GFelement& lhs, const GFelement& rhs)
{
    if (lhs.field != rhs.field) throw ErrorIncompatibleFields;

    return lhs.degree == rhs.degree;
}


bool operator==(Fint left, const GFelement &right) {
    return left == right.degree;
}


bool operator==(const GFelement &left, Fint right) {
    return left.degree == right;
}



//------------------------------------------------------
bool operator!=(const GFelement& lhs, const GFelement& rhs)
{
    return !(lhs == rhs);
}



bool operator!=(Fint left, const GFelement &right) {
    return !(left == right);
}


bool operator!=(const GFelement &left, Fint right) {
    return !(left == right);
}


//------------------------------------------------------
bool operator<(const GFelement& lhs, const GFelement& rhs)
{
    if (lhs.field != rhs.field) throw ErrorIncompatibleFields;

    return lhs.degree < rhs.degree;
}



//------------------------------------------------------
bool operator>(const GFelement& lhs, const GFelement& rhs)
{
   if (lhs.field != rhs.field) throw ErrorIncompatibleFields;

   return lhs.degree > rhs.degree;
}



//------------------------------------------------------
bool operator<=(const GFelement& lhs, const GFelement& rhs)
{
    return !(lhs > rhs);
}



//------------------------------------------------------
bool operator>=(const GFelement& lhs, const GFelement& rhs)
{
    return !(lhs < rhs);
}



//------------------------------------------------------
GFelement GFelement::operator+(const GFelement& rhs) const
{
    if (rhs.field != field) throw ErrorIncompatibleFields;

    GFelement result(rhs.field, rhs.field->add(degree, rhs.degree));
    return result;
}

GFelement GFelement::operator+(Fint rhs) {
    return *this + GFelement(this->field, rhs);
}

//------------------------------------------------------
GFelement GFelement::operator-(const GFelement& rhs) const
{
    if (rhs.field != field) throw ErrorIncompatibleFields;

    GFelement result(rhs.field, rhs.field->subtract(degree, rhs.degree));
    return result;
}

GFelement GFelement::operator-(Fint rhs) {
    return *this - GFelement(this->field, rhs);
}

//------------------------------------------------------
GFelement GFelement::operator*(const GFelement& rhs) const
{
//    auto if_this_gf_prime = typeid(this->field) == typeid(GaloisFieldPrime&);
//    auto if_rhs_gf_prime = typeid(rhs.field) == typeid(GaloisFieldPrime&);
//    auto if_this_gf_prime = dynamic_cast<GaloisFieldPrime*>(this->field); typeid(this->field) == typeid(GaloisFieldPrime&);
//    auto if_rhs_gf_prime = typeid(rhs.field) == typeid(GaloisFieldPrime&);
//
//    // if if_lhs_gf_prime (or if_rhs_gf_prime) is True, then GaloisFieldPrime, else GaloisField
//    if (!if_this_gf_prime && if_rhs_gf_prime) {
//        GFelement result(this->field, this->field->sum_times(this->degree, rhs.degree));
//        return result;
//    } else if (if_this_gf_prime && !if_rhs_gf_prime) {
//        GFelement result(rhs.field, rhs.field->sum_times(rhs.degree, this->degree));
//        return result;
//    }

    if (rhs.field != field) throw ErrorIncompatibleFields;

    GFelement result(rhs.field, rhs.field->multiply(degree, rhs.degree));
    return result;
}



//------------------------------------------------------
GFelement GFelement::operator/(const GFelement& rhs) const
{
    return  *this * rhs.inverse();
}



//------------------------------------------------------
GFelement GFelement::operator*(Fint rhs) const
{
    GFelement result(field, field->multiply(degree, rhs));
    return result;
}



//------------------------------------------------------
GFelement operator*(const GFelement& lhs, Fint rhs)
{
    GFelement result(lhs.field, lhs.field->multiply(lhs.degree, rhs));
    return result;
}



//------------------------------------------------------
GFelement operator/(const GFelement& lhs, Fint rhs)
{
    GFelement result(lhs.field, lhs.field->multiply(lhs.degree,
                                                   lhs.field->inverse(rhs)));
    return result;
}



//------------------------------------------------------
ostream& operator<<(ostream& output, const GFelement& rhs)
{
    output << rhs.degree;
    return output;
}



//------------------------------------------------------
GFelement GFelement::inverse() const
{
    GFelement result(field, field->inverse(degree));
    return result;
}



//------------------------------------------------------
GFelement GFelement::power(Fint p) const
{
    if (p < 0) throw std::logic_error("Power degree must be >= 0.");
    if (degree == 0 || degree == 1)
        return *this;

    Fint result_degree = 1;
    for (int i = 0; i < p; i++)
        result_degree = field->multiply(result_degree, degree);

    GFelement result(field, result_degree);
    return result;
}

//------------------------------------------------------
//GFelement GFelement::sum_times(Fint times) const
//{
//    if (times < 0) throw std::logic_error("Operand must be >= 0.");
//
//    Fint result_degree = 0;
//    for (int i = 1; i <= times; i++)
//        result_degree = field.add(result_degree, degree);
//
//    GFelement result(field, result_degree);
//    return result;
//}

GFelement GFelement::sum_times(Fint times) const {
    if (times < 0) throw std::logic_error("Operand must be >= 0.");

    auto t = times % field->getModulus();
    GFelement res(field, field->sum_times(degree, t));
    return res;
}

GFelement GFelement::summed_times(Fint times) const {
    if (times < 0) throw std::logic_error("Operand must be >= 0.");

    auto t = times % field->getModulus();
    GFelement res(field, field->summed_times(t, degree));
    return res;
}

GFelement GFelement::trace() const {
    GFelement tr(this->field, 0);
    for (auto i = 0; i < field->getDimension(); i++) {
        tr = tr + this->power(pow(this->field->getModulus(), i));
    }

    return tr;
}

/* TODO
 * 1) add fields compatibility check
 */
GFelement GFelement::dot(const std::vector<GFelement> &op1, const std::vector<GFelement> &op2) {
    if (op1.size() != op2.size())
        throw std::logic_error("Operands length must match");

    const GaloisField *field = op1[0].getField();

    GFelement s(field, 0);
    for (size_t i = 0; i < op1.size(); i++) {
        s = s + op1[i] * op2[i];
    }

    return s;
}

GFelement GFelement::dot(const std::vector<GFelement> &op1, const std::vector<int> &op2) {
    if (op1.size() != op2.size())
        throw std::logic_error("Operands length must match");

    auto field = op1[0].getField();

    GFelement s(field, 0);
    for (size_t i = 0; i < op1.size(); i++) {
//        std::cout << op1[i] << " * " << op2[i] << " = " << op1[i] * GFelement(field, op2[i]) << std::endl;
        s = s + op1[i] * GFelement(field, op2[i]);
    }

    return s;
}

GFelement GFelement::dot(const std::vector<int> &op1, const std::vector<GFelement> &op2) {
    return dot(op2, op1);
}

GFelement GFelement::dotint(const std::vector<int> &op1, const std::vector<GFelement> &op2) {
    if (op1.size() != op2.size())
        throw std::logic_error("Operands length must match");

    const GaloisField *field = op2[0].getField();

    GFelement s(field, 0);
    for (size_t i = 0; i < op1.size(); i++) {
//        std::cout << op1[i] << " ";
        s = s + op2[i].sum_times(op1[i]);
//        std::cout << s << " ";
    }
//    std::cout << std::endl;

    return s;
}

GFelement GFelement::dotint(const std::vector<GFelement> &op1, const std::vector<int> &op2) {
    return dotint(op2, op1);
}

std::vector<GFelement> GFelement::to_gf(const GaloisField *field, const vector<Fint> &vector) {
    std::vector<GFelement> res(vector.size());
    for (size_t i = 0; i < res.size(); i++) {
        res[i] = GFelement(field, vector[i]);
    }
    return res;
}

std::vector<Fint> GFelement::as_vector() const {
    return field->elems()[degree];
}

} // namespace

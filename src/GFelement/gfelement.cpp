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
using namespace std;

#include "gfelement.h"

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
    if (field != nullptr) {
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
GFelement GFelement::operator+(const GFelement& lhs) const
{
    if (lhs.field != field) throw ErrorIncompatibleFields;

    GFelement result(lhs.field, lhs.field->add(lhs.degree, degree));
    return result;
}



//------------------------------------------------------
GFelement GFelement::operator-(const GFelement& lhs) const
{
    if (lhs.field != field) throw ErrorIncompatibleFields;

    GFelement result(lhs.field, lhs.field->subtract(lhs.degree, degree));
    return result;
}



//------------------------------------------------------
GFelement GFelement::operator*(const GFelement& lhs) const
{
    if (lhs.field != field) throw ErrorIncompatibleFields;

    GFelement result(lhs.field, lhs.field->multiply(lhs.degree, degree));
    return result;
}



//------------------------------------------------------
GFelement GFelement::operator/(const GFelement& lhs) const
{
    return  lhs * this->inverse();
}



//------------------------------------------------------
GFelement GFelement::operator*(Fint lhs) const
{
    GFelement result(field, field->multiply(lhs, degree));
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
    if (degree == 0)
        return *this;

    Fint result_degree = 1;
    for (int i = 1; i < p; i++)
        result_degree = field->multiply(result_degree, degree);

    GFelement result(field, result_degree);
    return result;
}

//------------------------------------------------------
GFelement GFelement::sum_times(Fint times) const
{
    if (times < 0) throw std::logic_error("Operand must be >= 0.");

    Fint result_degree = 0;
    for (int i = 1; i <= times; i++)
        result_degree = field->add(result_degree, degree);

    GFelement result(field, result_degree);
    return result;
}


} // namespace

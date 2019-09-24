/*******************************************************************************
   Project: C++ Library for General Galois Field Arithmetic

   Language: C++ 2007	   
   Author: Saied H. Khayat
   Date:   Feb 2013
   URL: https://github.com/saiedhk
   
   Copyright Notice: Free use of this library is permitted under the
   guidelines and in accordance with the MIT License (MIT).
   http://opensource.org/licenses/MIT

*******************************************************************************/

#ifndef GFELEMENT_H
#define GFELEMENT_H

#include "GaloisField/galoisfield.h"
#include "Arith/polynomial_arith.h"

#include <iostream>

namespace galoiscpp
{

/**
   Defines a Galois Field element. Also defines arithmetic operations on Galois Field elements.

   In mathematics, an element of the Galois Field GF(p^k) is a polynomial of degree less than k
   with coefficients chosen from the set {0,1,2,..., p-1}.

   In this program, a GFelement contains the following private data:
   @param field a pointer to the Galois Field this GFelement belongs to
   @param modulus the prime modulus of the Galois Field
   @param dimension the dminesion of the Galois Field
   @param polynomial the array representing the polynomial describing this GFelement
*/
class GFelement
{

    public:
        GFelement() : degree(-1) {};
        GFelement(const GaloisField *gf, Fint degree) : field(gf), degree(degree) {};
        GFelement(const GFelement& gfe) : field(gfe.field), degree(gfe.degree) {};
        ~GFelement() {};

        GFelement& operator=(const GFelement& right);
        GFelement& operator=(Fint rhs);

        friend bool operator==(const GFelement& left, const GFelement& right);
        friend bool operator==(Fint left, const GFelement& right);
        friend bool operator==(const GFelement& left, Fint right);
        friend bool operator!=(const GFelement& left, const GFelement& right);
        friend bool operator!=(Fint left, const GFelement& right);
        friend bool operator!=(const GFelement& left, Fint right);

        friend bool operator<(const GFelement& left, const GFelement& right);
        friend bool operator>(const GFelement& left, const GFelement& right);
        friend bool operator<=(const GFelement& left, const GFelement& right);
        friend bool operator>=(const GFelement& left, const GFelement& right);

        GFelement operator+(const GFelement& rhs) const;
        GFelement operator+(Fint rhs);
        GFelement operator-(const GFelement& rhs) const;
        GFelement operator-(Fint rhs);
        GFelement operator-() const;

        GFelement operator*(const GFelement& rhs) const;
        GFelement operator*(Fint rhs) const;
        friend GFelement operator*(const GFelement& left, Fint right);
        GFelement operator/(const GFelement& rhs) const;
        friend GFelement operator/(const GFelement& left, Fint right);

        friend ostream& operator<<(ostream& output, const GFelement& right);

        GFelement inverse() const;

        GFelement power(Fint p) const;
        GFelement sum_times(Fint times) const;
//        GFelement sum_times_table(Fint times) const;
        GFelement summed_times(Fint times) const;
        GFelement trace() const;

        static GFelement dot(const std::vector<GFelement> &op1, const std::vector<GFelement> &op2);
        static GFelement dotint(const std::vector<int> &op1, const std::vector<GFelement> &op2);
        static GFelement dotint(const std::vector<GFelement> &op2, const std::vector<int> &op1);

        static std::vector<GFelement> to_gf(const GaloisField *field, const std::vector<Fint> &vector);
        std::vector<Fint> as_vector() const;

        const GaloisField* getField() const;
        Fint getDegree() const;

    private:
        const GaloisField *field;  // Galois Field associated with this polynomial
        Fint degree;
};



//------------------------------------------------
// Inline class functions
//------------------------------------------------

inline const GaloisField* GFelement::getField() const
{
   return field;
}


inline Fint GFelement::getDegree() const
{
   return degree;
}


} // namespace shk_galoiscpp

#endif

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

#include <iostream>
using namespace std;

#include "galoisfield.h"
#include "polynomial_arith.h"


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

        //---------------------------
        // constructors
        //---------------------------

        GFelement() : field(nullptr), degree(-1) {};

        /**
         Contructs GFelement without initializing its coefficients
         @param gf pointer to a Galois Field
        */
        GFelement(GaloisField* gf) : field(gf), degree(0) {};


        /**
         Contructs a GFelement and initializes its coefficients
         @param gf pointer to a Galois Field
         @param poly an external array containing coefficients (length of array = dimension)
        */
        GFelement(GaloisField* gf, Fint degree) : field(gf), degree(degree) {};


        /**
         Constructs and clones a GFelement.
         @param gfe reference to a GFelement
        */
        GFelement(const GFelement& gfe) : field(gfe.field), degree(gfe.degree) {};


        /**
         Destructs a GFelement and frees its allocated memory
        */
        ~GFelement() {field = nullptr; degree = -1;};


        //---------------------------
        // accessors
        //---------------------------

        /**
         Returns the pointer to the GF this GFelement is associated with.
        */
        GaloisField* getField() const;


        /**
         Returns the i-th coefficient of this GFelement's polynomial
         @param i the index to the i-th coefficient
        */
        Fint getDegree() const;


        //---------------------------
        // operators
        //---------------------------

        /**
         Assigns/copies the values of another GFelement to this GFelement
         @param right GFelement on the right of = sign
         @returns GFelement on the left of = sign (this GFelement)
        */
        GFelement& operator=(const GFelement& right);


        GFelement& operator=(Fint rhs);


        /**
         Returns the negated copy of this GFelement
        */
        GFelement operator-() const;


        /**
         Compares two GFelements, checks for equality
         @param right GFelement on the right of == sign
         @param left GFelement on the left of == sign
         @returns true if equal, false otherwise
        */
        friend bool operator==(const GFelement& left, const GFelement& right);

        friend bool operator==(Fint left, const GFelement& right);
        friend bool operator==(const GFelement& left, Fint right);

        /**
         Compares two GFelements, checks for inequality
         @param right GFelement on the right of != sign
         @param left GFelement on the left of != sign
         @returns true if equal, false otherwise
        */
        friend bool operator!=(const GFelement& left, const GFelement& right);
        friend bool operator!=(Fint left, const GFelement& right);
        friend bool operator!=(const GFelement& left, Fint right);


        /**
         Compares two GFelements based on their lexicographical ordering.
         @param right GFelement on the right of < sign
         @param left GFelement on the left of < sign
         @returns true if left < right, false otherwise
        */
        friend bool operator<(const GFelement& left, const GFelement& right);


        /**
         Compares two GFelements based on their lexicographical ordering.
         @param right GFelement on the right of > sign
         @param left GFelement on the left of > sign
         @returns true if left > right, false otherwise
        */
        friend bool operator>(const GFelement& left, const GFelement& right);


        /**
         Compares two GFelements based on their lexicographical ordering.
         @param right GFelement on the right of <= sign
         @param left GFelement on the left of <= sign
         @returns true if left <= right, false otherwise
        */
        friend bool operator<=(const GFelement& left, const GFelement& right);


        /**
         Compares two GFelements based on their lexicographical ordering.
         @param right GFelement on the right of >= sign
         @param left GFelement on the left of >= sign
         @returns true if left >= right, false otherwise
        */
        friend bool operator>=(const GFelement& left, const GFelement& right);


        /**
         Adds two GFelements
         @param right GFelement on the right of + sign
         @param left GFelement on the left of + sign
         @returns sum of two operands
        */
        GFelement operator+(const GFelement& lhs) const;


        /**
         Subtracts two GFelements
         @param right GFelement on the right of - sign
         @param left GFelement on the left of - sign
         @returns difference of two operands
        */
        GFelement operator-(const GFelement& lhs) const;


        /**
         Multiplies two GFelements
         @param right GFelement on the right of * sign
         @param left GFelement on the left of * sign
         @returns product of two operands
        */
        GFelement operator*(const GFelement& lhs) const;


        /**
         Divides two GFelements
         @param right GFelement on the right of / sign
         @param left GFelement on the left of / sign
         @returns quotient of two operands
        */
        GFelement operator/(const GFelement& lhs) const;


        /**
         Multiplies a GFelement by a scalar (i.e. a Fint)
         @param right GFelement on the right of * sign
         @param left Fint on the left of * sign
         @returns product of two operands
        */
        GFelement operator*(Fint left) const;


        /**
         Multiplies a GFelement by a scalar (i.e. a Fint)
         @param right Fint on the right of * sign
         @param left GFelement on the left of * sign
         @returns product of two operands
        */
        friend GFelement operator*(const GFelement& left, Fint right);


        /**
         Divides a GFelement by a scalar (i.e. a Fint)
         @param right Fint on the right of / sign
         @param left GFelement on the left of / sign
         @returns quotient of two operands
         @throws ErrorDivideByZero
        */
        friend GFelement operator/(const GFelement& left, Fint right);


        /**
         Outputs a GFelement to the standard output
         @param right GFelement on the right of << sign
         @param output an output stream on the left of << sign
         @returns output stream
        */
        friend ostream& operator<<(ostream& output, const GFelement& right);


        /**
         Determines the multiplicative inverse of a GFelement.
         @param gfe input GFelement
         @returns GFelement (the multiplicative inverse of gfe)
         @throws ErrorNoInverse, ErrorInvalidArraySize
        */
        GFelement inverse() const;


        /**
         Raises a GFelement G to a positive integer power m

         It uses the repeated square-and-multiply algorithm for exponentiation in GF(p^k).

         (See Chapter 2, Algorithm 2.227 in Handbook of Applied Cryptography,
         by Menezes, van Oorschot, Vanstone, CRC Press, 1996.)

         @param G an input GFelement
         @param m an input positive integer
         @returns G raised to the power m
         @throws ErrorInvalidArraySize, ErrorDivideByZero, ErrorIncompatibleFields
        */
        GFelement power(Fint p) const;

        GFelement sum_times(Fint times) const;


    private:
        GaloisField* field;  // Galois Field associated with this polynomial
        Fint degree;
};



//------------------------------------------------
// Inline class functions
//------------------------------------------------

inline GaloisField* GFelement::getField() const
{
   return field;
}


inline Fint GFelement::getDegree() const
{
   return degree;
}


} // namespace shk_galoiscpp

#endif

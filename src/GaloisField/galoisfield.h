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

#ifndef GALOISFIELD_H
#define GALOISFIELD_H

#include "typedefs.h"
#include <cassert>
#include <vector>



namespace galoiscpp
{

/**
   This class defines a Galois Field of characteristic p and dimesion k (for GF(p^k)).
   It stores the following (private) data about the field:
   @param modulus The modulus (aka the characteristics p), which MUST be a prime number.
   @param dimension The dimension (k) of Galois Field (a nonzero positive integer).
   @param reductpoly A polynomial of degree k over the prime field GF(p). This polynomial
   MUST be irreducible. The program does not check this. YOU MUST correctly choose
   this polynomial.

   The attributes of a GaloisField are set when it is instantiated.
   Once instantiated, the attributes of cannot be changed.
*/
class GaloisField
{
   public:

      //---------------------------
      // constructors
      //---------------------------
      /**
         Constructs a Galois Field (GF)
         @param mod a positive prime integer (the modulus of GF)
         @param dim a positive integer (the dimension of GF)
         @param poly an array of length dim(the reduction polynomial of GF)
         (must be an irreducible polynomial over the prime field of GF)
      */
      GaloisField(Fint mod, Int dim, Fint poly[]);

      GaloisField(Fint mod, Int dim);

      /**
         Destructs GaloisField, freeing allocated space for reduction polynomial in it.
      */
      ~GaloisField();

      //---------------------------
      // accessors
      //---------------------------

      /**
         Returns modulus of GF.
      */
      Fint getModulus() const;

      /**
         Returns dimension of GF
      */
      Int getDimension() const;

      /**
         Returns i-th coefficient of reduction polynomial of GF
      */
      Fint reductionPolynomial(Int i) const;

      /**
         Outputs information about this Galois Field to the standard output.
         @param gf GaloisField on the right of << sign
         @param output an output stream on the left of << sign
         @returns output stream
      */
      friend std::ostream& operator<<(std::ostream& output, const GaloisField& gf);

      // Getters
      std::vector<std::vector<Fint>> elems();
      Fint add(Fint lhs, Fint rhs);
      Fint subtract(Fint lhs, Fint rhs);
      Fint multiply(Fint lhs, Fint rhs);
      Fint inverse(Fint op);

      Fint get_size();

      /* TODO
       * 1) get rid of all C-like stuff
       */
   private:
      Fint modulus;      // prime modulus (characteristic) of Galois Field
      Int dimension;     // dimension of Galois Field
      Fint* reductpoly;  // reduction polynomial of Galois Field
      Fint size;

      std::vector<std::vector<Fint>> field_table;    // correspondence between a power of a primitive elem and a poly
      std::vector<std::vector<Fint>> add_table;             // a table for addition
      std::vector<std::vector<Fint>> sub_table;             // a table for subtraction
      std::vector<std::vector<Fint>> mult_table;            // a table for multiplication
      std::vector<Fint> inv_table;             // table inverse elems

      void create_tables();
};



//--------------------------------------
// Inline class functions
//--------------------------------------

inline Fint GaloisField::getModulus() const
{
   return modulus;
}


inline Int GaloisField::getDimension() const
{
   return dimension;
}


inline Fint GaloisField::reductionPolynomial(Int i) const
{
   return reductpoly[i];
}


} // namespace shk_galoiscpp

#endif

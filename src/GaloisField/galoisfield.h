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

#include <ostream>
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
      GaloisField() : modulus(-1), dimension(-1), reductpoly(nullptr), size(-1) {};

      /**
         Constructs a Galois Field (GF)
         @param mod a positive prime integer (the modulus of GF)
         @param dim a positive integer (the dimension of GF)
         @param poly an array of length dim(the reduction polynomial of GF)
         (must be an irreducible polynomial over the prime field of GF)
      */
      GaloisField(Fint mod, Int dim, Fint poly[]);

      GaloisField(Fint mod, Int dim);

      GaloisField(const GaloisField &gf);

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
      virtual Fint getModulus() const;

      /**
         Returns dimension of GF
      */
      virtual Int getDimension() const;

      /**
         Returns i-th coefficient of reduction polynomial of GF
      */
      Fint reductionPolynomial(Int i) const;

      Fint get_generator() const { return generator; }

      /**
         Outputs information about this Galois Field to the standard output.
         @param gf GaloisField on the right of << sign
         @param output an output stream on the left of << sign
         @returns output stream
      */
      friend std::ostream& operator<<(std::ostream& output, const GaloisField& gf);

      friend bool operator==(const GaloisField &lhs, const GaloisField &rhs);
      friend bool operator!=(const GaloisField &lhs, const GaloisField &rhs);

      GaloisField& operator=(const GaloisField &rhs);

      // Getters
      std::vector<std::vector<Fint>> elems() const;
      Fint add(Fint lhs, Fint rhs) const;
      Fint subtract(Fint lhs, Fint rhs) const;
      Fint multiply(Fint lhs, Fint rhs) const;
      Fint inverse(Fint op) const;

      Fint sum_times(Fint op, int times) const;
      Fint summed_times(Fint times, Fint res) const;

      std::vector<std::vector<Fint>> *field_table_ptr() { return &field_table; } // TODO: !

      virtual Fint get_size() const;

      std::vector<Fint> find_subfield(int sub_m) const;

      /* TODO
       * 1) get rid of all C-like stuff
       */
    private:
      std::vector<std::vector<Fint>> field_table;    // correspondence between a power of a primitive elem and a poly

    protected:
      Fint modulus;      // prime modulus (characteristic) of Galois Field
      Int dimension;     // dimension of Galois Field
      Fint* reductpoly;  // reduction polynomial of Galois Field
      Fint size;

      Fint generator;

      std::vector<std::vector<Fint>> add_table;             // a table for addition
      std::vector<std::vector<Fint>> sub_table;             // a table for subtraction
      std::vector<std::vector<Fint>> mult_table;            // a table for multiplication
      std::vector<Fint> inv_table;             // table inverse elems

      std::vector<std::vector<Fint>> sum_times_table;
      std::vector<std::vector<Fint>> sum_times_table_transposed;

      virtual void create_tables();
};

inline Fint GaloisField::getModulus() const { return modulus; }
inline Int GaloisField::getDimension() const { return dimension; }
inline Fint GaloisField::reductionPolynomial(Int i) const { return reductpoly[i]; }

} // namespace shk_galoiscpp

#endif

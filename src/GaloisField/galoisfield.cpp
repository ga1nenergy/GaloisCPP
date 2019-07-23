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
#include <vector>
#include <cmath>
#include <cstring>
using namespace std;

#include <cassert>
#include "galoisfield.h"
#include "polynomial_arith.h"
#include "irreducible_poly.hpp"

namespace galoiscpp
{

//------------------------------------------------------
GaloisField::GaloisField(Fint mod, Int dim, Fint poly[])
{
   assert (mod>1);
   assert (dim>0);

   modulus = mod;
   dimension = dim;
   reductpoly = new Fint[dimension+1];
   size = Fint(pow(modulus, dimension)) - 1;

   for (Int i=0; i < dimension+1; i++)
   {
      reductpoly[i] = poly[i];
   }

   create_tables();
}


//------------------------------------------------------
GaloisField::GaloisField(Fint mod, Int dim) {
   assert (mod>1);
   assert (dim>0);

   modulus = mod;
   dimension = dim;
   reductpoly = new Fint[dimension+1];
   size = Fint(pow(modulus, dimension));

   get_irreducible_poly(modulus, dimension, reductpoly);
   create_tables();
}


//------------------------------------------------------
GaloisField::~GaloisField()
{
   delete [] reductpoly;
}


//------------------------------------------------------
ostream& operator<<(ostream& output, const GaloisField& gf)
{
   output << "\nPrime Modulus: " << gf.modulus;
   output << "\nDimension: " << gf.dimension;

   output << "\nReduction Polynomial Coefficients: (";
   for (Int i=0; i<=gf.dimension; i++)
   {
      if (i < gf.dimension)
         output << gf.reductpoly[i] << " ";
      else
         output << gf.reductpoly[i] << ")\n";
   }

    output << "Field:" << endl;
    for (int i = 0; i < gf.size; i++) {
        output << i << ": ";
        for (int j = 0; j < gf.dimension; j++)
            output << gf.field_table[i][j];
        output << endl;
    }

    output << endl << "Addition table" << endl;
    for (auto & v : gf.add_table) {
        for (auto &e : v)
            output << e << " ";
        output << endl;
    }

    output << endl << "Subtraction table" << endl;
    for (auto & v : gf.sub_table) {
        for (auto &e : v)
            output << e << " ";
        output << endl;
    }

    output << endl << "Multiplication table: " << endl;
    for (auto & v : gf.mult_table) {
        for (auto &e : v)
            output << e << " ";
        output << endl;
    }

    output << endl << "Inverse table: " << endl;
    for (int i = 0; i < gf.size; i++) {
        output << i << ": " << gf.inv_table[i] << endl;
    }

    return output;
}

void GaloisField::create_tables() {
    // 1) create field
    field_table = std::vector<std::vector<Fint>>(size);
    field_table[0] = std::vector<Fint>(dimension, 0);
    for (int i = 1; i < size; i++) {
        Fint* poly = new Fint[i]; memset(poly, 0, i * sizeof(Fint)); poly[i - 1] = 1;
        field_table[i] = std::vector<Fint>(dimension, 0);
        Fint* Q = new Fint[i];
        Fint* R = new Fint[i];

        polynomialDivide(poly, i, reductpoly, dimension + 1, Q, i, R, i, modulus);

        memcpy(field_table[i].data(), R, dimension * sizeof(Fint)); // sic!
        delete [] poly; delete [] Q; delete [] R;
    }

    // 2) create addition table
    add_table = std::vector<std::vector<Fint>>(size);
    for (int i = 0; i < size; i++) {
        add_table[i] = std::vector<Fint>(size);
        for (int j = 0 ; j < size; j++) {
            Fint* res = new Fint[dimension];
            polynomialAdd(field_table[i].data(), dimension,
                       field_table[j].data(), dimension,
                       res, dimension, modulus);

            for (int k = 0; k < size; k++) {
                if (!memcmp(res, field_table[k].data(), dimension * sizeof(Fint))) {
                    add_table[i][j] = k;
                    break;
                }
            }

            delete [] res;
        }
    }

    // 3) create subtraction table
    sub_table = std::vector<std::vector<Fint>>(size);
    for (int i = 0; i < size; i++) {
        sub_table[i] = std::vector<Fint>(size);
        for (int j = 0; j < size; j++) {
            Fint *res = new Fint[dimension];
            polynomialSubtract(field_table[i].data(), dimension,
                               field_table[j].data(), dimension,
                               res, dimension, modulus);

            for (int k = 0; k < size; k++) {
                if (!memcmp(res, field_table[k].data(), dimension * sizeof(Fint))) {
                    sub_table[i][j] = k;
                    break;
                }
            }

            delete[] res;
        }
    }

    // 3) create multiplication table
    mult_table = std::vector<std::vector<Fint>>(size);
    mult_table[0] = std::vector<Fint>(size, 0);
    for (int i = 1; i < size; i++) {
        mult_table[i] = std::vector<Fint>(size, 0);
        for (int j = 1; j < size; j++) {
            mult_table[i][j] = (i + j - 2) % (size - 1) + 1;
        }
    }

    // check
//    cout << std::endl << "Multiplication table check: " << std::endl;
//    for (int i = 1; i < size; i++) {
//        for (int j = 1; j < size; j++) {
//            Fint *res = new Fint[dimension*2];
//            polynomialMultiply(field_table[i].data(), dimension,
//                               field_table[j].data(), dimension,
//                               res, 2*dimension, modulus);
//
//            Fint *res_rem = new Fint[dimension + 1];
//            Fint* Q = new Fint[dimension + 1];
//
//            polynomialDivide(res, dimension + 1, reductpoly, dimension + 1, Q, dimension + 1, res_rem, dimension + 1, modulus);
//
//            for (int k = 0; k < size; k++) {
//                if (!memcmp(res_rem, field_table[k].data(), dimension * sizeof(Fint))) {
//                    std::cout << k << " ";
//                }
//            }
//
//            delete[] res; delete[] res_rem; delete[] Q;
//        }
//        std::cout << std::endl;
//    }

    // 4) create inverse table
    inv_table = std::vector<Fint>(size, 0);

    for (int i = 1; i < size; i++) {
        inv_table[i] = (size - i) % (size - 1) + 1;
    }
}

Fint GaloisField::add(Fint lhs, Fint rhs) {
    return add_table[lhs][rhs];
}

Fint GaloisField::subtract(Fint lhs, Fint rhs) {
    return sub_table[lhs][rhs];
}

Fint GaloisField::multiply(Fint lhs, Fint rhs) {
    return mult_table[lhs][rhs];
}

Fint GaloisField::inverse(Fint op) {
    return inv_table[op];
}

std::vector<std::vector<Fint>> GaloisField::elems() {
    return field_table;
}

Fint GaloisField::get_size() {
    return size;
}

} // namespace shk_galoiscpp

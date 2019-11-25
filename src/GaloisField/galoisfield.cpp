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
#include <algorithm>
#include <iterator>
using namespace std;

#include <cassert>
#include "galoisfield.h"
#include "Arith/polynomial_arith.h"
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

   std::cout << "Field is created!. Address: " << this << std::endl;
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
            output << gf.field_table[i][j] << " ";
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

    output << endl << "Sum times table: " << endl;
    for (auto const &row : gf.sum_times_table) {
        std::copy(row.begin(), row.end(), std::ostream_iterator<int>(std::cout, " "));
        std::cout << std::endl;
    }

    output << endl << "Summed times table: " << endl;
    for (auto const &row : gf.sum_times_table_transposed) {
        std::copy(row.begin(), row.end(), std::ostream_iterator<int>(std::cout, " "));
        std::cout << std::endl;
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

        for (int j = 0; j < i && j < dimension; j++) {
            field_table[i][j] = R[j];
        }
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

    // 4) create inverse table
    inv_table = std::vector<Fint>(size, 0);

    for (int i = 1; i < size; i++) {
        inv_table[i] = (size - i) % (size - 1) + 1;
    }

    int max = modulus;
    for (auto i = 0; i < size; i++) {
        std::vector<Fint> row(max, 0);
        for (auto j = 0; j < max; j++) {
            for (auto k = 1; k <= j; k++) {
                row[j] = add(row[j], i);
            }
        }
        sum_times_table.push_back(row);
    }

    for (auto i = 0; i < max; i++) {
        std::vector<Fint> row(size);
        for (auto j = 0; j < size; j++) {
            row[j] = sum_times_table[j][i];
        }
        sum_times_table_transposed.push_back(row);
    }
}

Fint GaloisField::add(Fint lhs, Fint rhs) const {
    return add_table[lhs][rhs];
}

Fint GaloisField::subtract(Fint lhs, Fint rhs) const {
    return sub_table[lhs][rhs];
}

Fint GaloisField::multiply(Fint lhs, Fint rhs) const {
    return mult_table[lhs][rhs];
}

Fint GaloisField::inverse(Fint op) const {
    return inv_table[op];
}

std::vector<std::vector<Fint>> GaloisField::elems() const {
    return field_table;
}

Fint GaloisField::get_size() const {
    return size;
}

/* TODO
 * 1) think how to combine instances of a field and its subfield
 * 2) improve "while-loop" logic
 */
std::vector<Fint> GaloisField::find_subfield(int sub_m) const {
    if (sub_m < 1) {
        throw std::logic_error("Field extension must be a natural number");
    }

    Fint subfield_size = pow(modulus, sub_m) - 1;
    std::vector<Fint> set;
    for (Fint i = 2; i < size; i++) {
        set.clear();
        Fint elem = 1;

        do {
            elem = multiply(elem, i);
            set.push_back(elem);
        } while (set.size() != subfield_size); //!
        if (set[set.size() - 1] == 1)
            break;
    }

    return set;
}

    bool operator==(const GaloisField &lhs, const GaloisField &rhs) {
        return lhs.modulus == rhs.modulus && lhs.dimension == rhs.dimension && lhs.size == rhs.size &&
               !memcmp(lhs.reductpoly, rhs.reductpoly, (lhs.dimension + 1) * sizeof(Fint));
    }

    bool operator!=(const GaloisField &lhs, const GaloisField &rhs) {
        return !(lhs == rhs);
    }

    GaloisField& GaloisField::operator=(const GaloisField &gf) {
        modulus = gf.modulus;
        dimension = gf.dimension;

        reductpoly = new Fint[dimension + 1];
        memcpy(reductpoly, gf.reductpoly, (dimension + 1) * sizeof(Fint));

        size = gf.size;
        field_table = gf.field_table;
        add_table = gf.add_table;
        sub_table = gf.sub_table;
        mult_table = gf.mult_table;
        inv_table = gf.inv_table;
        sum_times_table = gf.sum_times_table;
        sum_times_table_transposed = gf.sum_times_table_transposed;

        return *this;
    }

    GaloisField::GaloisField(const GaloisField &gf) {
        modulus = gf.modulus;
        dimension = gf.dimension;

        reductpoly = new Fint[dimension + 1];
        memcpy(reductpoly, gf.reductpoly, (dimension + 1) * sizeof(Fint));

        size = gf.size;
        field_table = gf.field_table;
        add_table = gf.add_table;
        sub_table = gf.sub_table;
        mult_table = gf.mult_table;
        inv_table = gf.inv_table;
        sum_times_table = gf.sum_times_table;
        sum_times_table_transposed = gf.sum_times_table_transposed;
    }

    Fint GaloisField::sum_times(Fint op, int times) const {
        return sum_times_table[op][times];
    }

    Fint GaloisField::summed_times(Fint times, Fint res) const {
        auto elem = std::find(sum_times_table_transposed[times].begin(),
                              sum_times_table_transposed[times].end(),
                              res);
        if (elem != sum_times_table_transposed[times].end()) {
            return (elem - sum_times_table_transposed[times].begin());
        } else {
            throw std::logic_error("elem not found");
        }
    }

} // namespace shk_galoiscpp

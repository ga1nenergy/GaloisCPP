//
// Created by ga1nenergy on 11.11.2019.
//

#include <GaloisCPP/Arith/polynomial_arith.h>
#include "galois_field_prime.h"
#include <numeric>
#include <iterator>

namespace galoiscpp {
    GaloisFieldPrime::GaloisFieldPrime(Fint mod) {
        modulus = mod;
        size = mod;
        dimension = 1;
        reductpoly = new Fint[2]; // it is not used. added for compatibility with GaloisField

        create_tables();

        for (auto i = 2; i < size; i++) {
            auto counter = 1;
            auto buf = i;
            while (buf != 1) {
                buf = mult_table[buf][i];
                counter++;
            }

            if (counter == size - 1) {
                generator = i;
                break;
            }
        }

        std::cout << "Prime field is created!. Address: " << this << std::endl;
    }

    void GaloisFieldPrime::create_tables() {
        // 1) create field
        field_table = std::vector<Fint>(size);
        std::iota(field_table.begin(), field_table.end(), 0);

        // 2) create addition table
        add_table = std::vector<std::vector<Fint>>(size);
        for (int i = 0; i < size; i++) {
            add_table[i] = std::vector<Fint>(size);
            for (int j = 0 ; j < size; j++) {
                add_table[i][j] = (field_table[i] + field_table[j]) % modulus;
            }
        }

        // 3) create subtraction table
        sub_table = std::vector<std::vector<Fint>>(size);
        for (int i = 0; i < size; i++) {
            sub_table[i] = std::vector<Fint>(size);
            for (int j = 0; j < size; j++) {
                sub_table[i][j] = (field_table[i] + modulus - field_table[j]) % modulus;
            }
        }

        // 3) create multiplication table
        mult_table = std::vector<std::vector<Fint>>(size);
        mult_table[0] = std::vector<Fint>(size, 0);
        for (int i = 1; i < size; i++) {
            mult_table[i] = std::vector<Fint>(size, 0);
            for (int j = 1; j < size; j++) {
                mult_table[i][j] = (field_table[i] * field_table[j]) % modulus;
            }
        }

        // 4) create inverse table
        inv_table = std::vector<Fint>(size, 0);

        for (int i = 1; i < size; i++) {
            for (int j = 1; j < size; j++) {
                if (mult_table[i][j] == 1) {
                    inv_table[i] = j;
                    break;
                }
            }
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

    GaloisFieldPrime::GaloisFieldPrime(const GaloisFieldPrime &gf) {
        modulus = gf.modulus;

        size = gf.size;
        field_table = gf.field_table;
        add_table = gf.add_table;
        sub_table = gf.sub_table;
        mult_table = gf.mult_table;
        inv_table = gf.inv_table;
        sum_times_table = gf.sum_times_table;
        sum_times_table_transposed = gf.sum_times_table_transposed;
    }

    GaloisFieldPrime &GaloisFieldPrime::operator=(const GaloisFieldPrime &gf) {
        modulus = gf.modulus;

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

    std::ostream &operator<<(std::ostream &output, const GaloisFieldPrime &gf) {
        output << "\nCharacteristic: " << gf.modulus;

        output << endl << "Field: ";
        for (int i = 0; i < gf.size; i++) {
            output << i << " ";
        }

        output << endl << endl << "Addition table" << endl;
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


}
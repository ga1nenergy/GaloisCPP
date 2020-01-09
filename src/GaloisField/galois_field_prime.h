//
// Created by ga1nenergy on 11.11.2019.
//

#ifndef REMAINDER_CODE_LIST_DECODING_GALOIS_FIELD_PRIME_H
#define REMAINDER_CODE_LIST_DECODING_GALOIS_FIELD_PRIME_H

#include "galoisfield.h"
#include "typedefs.h"

namespace galoiscpp {
    class GaloisFieldPrime : public GaloisField {
    private:
        std::vector<Fint> field_table;

        void create_tables() override;

    public:
        GaloisFieldPrime(Fint mod);
        GaloisFieldPrime(const GaloisFieldPrime &field);
//        ~GaloisFieldPrime() {};

        friend std::ostream& operator<<(std::ostream& output, const GaloisFieldPrime& gf);

        friend bool operator==(const GaloisFieldPrime &lhs, const GaloisFieldPrime &rhs);
        friend bool operator!=(const GaloisFieldPrime &lhs, const GaloisFieldPrime &rhs);

        GaloisFieldPrime& operator=(const GaloisFieldPrime &rhs);

//        Fint getModulus() const override { return modulus; }
//        Fint get_size() const override { return size; }
    };
}


#endif //REMAINDER_CODE_LIST_DECODING_GALOIS_FIELD_PRIME_H

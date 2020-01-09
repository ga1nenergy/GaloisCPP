//
// Created by ga1nenergy on 03.01.2020.
//

#ifndef BCH_STRUCTURES_H
#define BCH_STRUCTURES_H

#include <cstddef>

template<typename T>
struct array_t {
    T *array;
    size_t size;

    T &operator[](int idx) {
        return array[idx];
    }
};

#endif //BCH_STRUCTURES_H

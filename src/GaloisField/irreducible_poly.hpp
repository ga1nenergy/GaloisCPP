#include "typedefs.h"
#include <cstring>
#include <vector>

namespace galoiscpp {
    void get_irreducible_poly(Fint mod, Int dim, Fint *pointer) {
        std::vector<Fint> poly;
        switch (mod) {
            case 2:
                switch (dim) {
                    case 1:
                        poly = {1, 1};
                        break;
                    case 2:
                        poly = {1, 1, 1};
                        break;
                    case 3:
                        poly = {1, 1, 0, 1};
                        break;
                    case 4:
                        poly = {1, 1, 0, 0, 1};
                        break;
                    case 5:
                        poly = {1, 0, 1, 0, 0, 1};
                        break;
                    case 6:
                        poly = {1, 1, 0, 0, 0, 0, 1};
                        break;
                    case 7:
                        poly = {1, 0, 0, 1, 0, 0, 0, 1};
                        break;
                    case 8:
                        poly = {1, 0, 1, 1, 1, 0, 0, 0, 1};
                        break;
                    case 9:
                        poly = {1, 0, 0, 0, 1, 0, 0, 0, 0, 1};
                        break;
                    case 10:
                        poly = {1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};
                }
                break;
            case 3:
                switch (dim) {
                    case 1:
                        poly = {1, 1};
                        break;
                    case 2:
                        poly = {1, 1, 2};
                        break;
                }
                break;
            case 5:
                switch (dim) {
                    case 1:
                        poly = {1, 2};
                        break;
                    case 2:
                        poly = {1, 1, 2};
                        break;
                }
                break;
            case 7:
                switch (dim) {
                    case 1:
                        poly = {1, 2};
                        break;
                    case 2:
                        poly = {1, 1, 3};
                        break;
                }
                break;
                //...
        }
        memcpy(pointer, poly.data(), (dim + 1) * sizeof(Fint));
    }
}
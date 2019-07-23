#include <iostream>
#include <vector>
#include <GFpoly/representations.h>
#include "gfelement.h"
#include "galoisfield.h"
#include "gfpoly.h"

using namespace galoiscpp;

int main() {

    GaloisField field(2, 3);
    std::cout << field << std::endl;
    std::vector<Fint> coefs = {4, 2, 1, 2, 3, 4};
    GFpoly op1(&field, coefs);
    coefs = {5, 6, 7};
    GFpoly op2(&field, coefs);

    coefs = {0, 0, 0, 0, 1};
    GFpoly op3(&field, coefs);

    std::cout << "op1: " << op1 << std::endl;
    std::cout << "op2: " << op2 << std::endl;

    auto inv_op2 = op2.inverse(op3);
    std::cout << "inverse: " << inv_op2 << std::endl;
    cout << "op2 * inv_op2: " << op2 * inv_op2 << std::endl;
    auto qr = (op2 * inv_op2) / op3;
    std::cout << "r: " << qr[1] << std::endl;

    GFpoly res = op1 + op1;
    qr = op1 / op2;
    auto prod = op1 * op2;
    std::cout << "prod: " << prod << std::endl;
    std::cout << qr[0] << std::endl;
    std::cout << qr[1] << std::endl;

    std::vector<Fint> v = {0, 1, 2, 3, 4, 5};
    GFpoly locator_poly = roots_to_poly(&field, v);
    std::cout << "locator_poly: " << locator_poly << std::endl;
    auto rts = locator_poly.roots();
    for (auto & rt : rts) {
        std::cout << rt << " ";
    }
    std::cout << std::endl;

    std::cout << "Hello, World!" << std::endl;
    return 0;
}
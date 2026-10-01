#include <iostream>
#include "ns2d/ns2d.h"

int main() {
    std::cout << "Hello, World!" << std::endl;

    if (checkSrc2d() == 0) {
        std::cout << "checkSrc2d() returned 0" << std::endl;
    } else {
        std::cout << "checkSrc2d() did not return 0" << std::endl;
    }

    return 0;
}
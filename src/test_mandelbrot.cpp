#include <iostream>
#include <cassert>
#include "rff2/calc/fp_complex_calculator.h"
#include "rff2/calc/fp_decimal_calculator.h"

using namespace merutilm::rff2;

void test_square() {
    // Test 1: (1 + 1i)^2 = 2i
    fp_complex_calculator z1("1", "1", 100);
    z1.square();
    std::cout << "z1^2 real: " << z1.getReal().double_value() << std::endl;
    std::cout << "z1^2 imag: " << z1.getImag().double_value() << std::endl;
    
    // Allow small float error due to double conversion for printing, but logic is integer based mostly
    assert(std::abs(z1.getReal().double_value() - 0.0) < 1e-10);
    assert(std::abs(z1.getImag().double_value() - 2.0) < 1e-10);

    // Test 2: (2 + 3i)^2 = -5 + 12i
    fp_complex_calculator z2("2", "3", 100);
    z2.square();
    std::cout << "z2^2 real: " << z2.getReal().double_value() << std::endl;
    std::cout << "z2^2 imag: " << z2.getImag().double_value() << std::endl;

    assert(std::abs(z2.getReal().double_value() - (-5.0)) < 1e-10);
    assert(std::abs(z2.getImag().double_value() - 12.0) < 1e-10);

    std::cout << "All tests passed!" << std::endl;
}

int main() {
    test_square();
    return 0;
}

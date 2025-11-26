//
// Created by Merutilm on 2025-05-05.
// Created by Super Fractal on 2025-11-25.
//

#pragma once

#include "dex.h"
#include <gmp.h>
#include <string>

namespace merutilm::rff2 {
    struct fp_decimal_calculator final {
        int exp2;
        mpz_t value;

        // Constructors
        fp_decimal_calculator() noexcept;
        explicit fp_decimal_calculator(mpz_srcptr value, int exp2);
        explicit fp_decimal_calculator(const dex &d, int exp10);
        explicit fp_decimal_calculator(double d, int exp10);
        explicit fp_decimal_calculator(const std::string &str, int exp10);
        ~fp_decimal_calculator();

        // Copy/Move semantics
        fp_decimal_calculator(const fp_decimal_calculator &other);
        fp_decimal_calculator &operator=(const fp_decimal_calculator &other);
        fp_decimal_calculator(fp_decimal_calculator &&other) noexcept;
        fp_decimal_calculator &operator=(fp_decimal_calculator &&other) noexcept;

        // Static arithmetic operations (assumes same exponents)
        static void fp_mul(fp_decimal_calculator &out, const fp_decimal_calculator &a, const fp_decimal_calculator &b) noexcept;
        static void fp_add(fp_decimal_calculator &out, const fp_decimal_calculator &a, const fp_decimal_calculator &b) noexcept;
        static void fp_sub(fp_decimal_calculator &out, const fp_decimal_calculator &a, const fp_decimal_calculator &b) noexcept;
        static void fp_div(fp_decimal_calculator &out, const fp_decimal_calculator &a, const fp_decimal_calculator &b) noexcept;
        static void fp_dbl(fp_decimal_calculator &out, const fp_decimal_calculator &target) noexcept;
        static void fp_hlv(fp_decimal_calculator &out, const fp_decimal_calculator &target) noexcept;
        static void fp_swap(fp_decimal_calculator &a, fp_decimal_calculator &b) noexcept;
        static void negate(fp_decimal_calculator &target) noexcept;

        // Utility functions
        static int exp2ToExp10(int exp2) noexcept;
        static int exp10ToExp2(int precision) noexcept;

        // Conversion
        [[nodiscard]] double double_value() const noexcept;
        void double_exp_value(dex *result) const;
        void setExp10(int exp10);

    private:
        // Thread-local temp storage for intermediate calculations
        static mpz_t& get_temp() noexcept;
    };
}
//
// Created by Merutilm on 2025-05-05.
// Created by Super Fractal on 2025-11-25.
//

#include "fp_decimal_calculator.h"
#include <cmath>
#include <cstring>
#include <limits>

// C++17 compatible bit_cast replacement
namespace {
    template<typename To, typename From>
    inline To bit_cast_compat(const From& src) noexcept {
        static_assert(sizeof(To) == sizeof(From), "Size mismatch");
        To dst;
        std::memcpy(&dst, &src, sizeof(To));
        return dst;
    }

    // Constants
    constexpr double LOG10_2 = 0.30102999566398119521373889472449;
}

namespace merutilm::rff2 {

    // Thread-local temp storage to avoid per-instance allocation
    mpz_t& fp_decimal_calculator::get_temp() noexcept {
        thread_local struct TempHolder {
            mpz_t t;
            TempHolder() { mpz_init(t); }
            ~TempHolder() { mpz_clear(t); }
        } holder;
        return holder.t;
    }

    fp_decimal_calculator::fp_decimal_calculator() noexcept : exp2(0) {
        mpz_init(value);
    }

    fp_decimal_calculator::fp_decimal_calculator(const mpz_srcptr val, const int e2) : exp2(e2) {
        mpz_init_set(value, val);
    }

    fp_decimal_calculator::fp_decimal_calculator(const dex &d, const int exp10) {
        mpz_init(value);
        exp2 = exp10ToExp2(exp10);

        mpf_t d1;
        mpf_init2(d1, static_cast<mp_bitcnt_t>(-exp2));
        mpf_set_d(d1, d.get_mantissa());

        const int mul_div = exp2 - d.get_exp2();
        if (mul_div < 0) {
            mpf_mul_2exp(d1, d1, static_cast<mp_bitcnt_t>(-mul_div));
        } else if (mul_div > 0) {
            mpf_div_2exp(d1, d1, static_cast<mp_bitcnt_t>(mul_div));
        }

        mpz_set_f(value, d1);
        mpf_clear(d1);
    }

    fp_decimal_calculator::fp_decimal_calculator(const double d, const int exp10) {
        mpz_init(value);
        exp2 = exp10ToExp2(exp10);

        mpf_t d1;
        mpf_init2(d1, static_cast<mp_bitcnt_t>(-exp2));
        mpf_set_d(d1, d);

        if (exp2 < 0) {
            mpf_mul_2exp(d1, d1, static_cast<mp_bitcnt_t>(-exp2));
        } else if (exp2 > 0) {
            mpf_div_2exp(d1, d1, static_cast<mp_bitcnt_t>(exp2));
        }

        mpz_set_f(value, d1);
        mpf_clear(d1);
    }

    fp_decimal_calculator::fp_decimal_calculator(const std::string &str, const int exp10) {
        mpz_init(value);
        exp2 = exp10ToExp2(exp10);

        mpf_t d1;
        mpf_init2(d1, static_cast<mp_bitcnt_t>(-exp2));
        mpf_set_str(d1, str.c_str(), 10);

        if (exp2 < 0) {
            mpf_mul_2exp(d1, d1, static_cast<mp_bitcnt_t>(-exp2));
        } else if (exp2 > 0) {
            mpf_div_2exp(d1, d1, static_cast<mp_bitcnt_t>(exp2));
        }

        mpz_set_f(value, d1);
        mpf_clear(d1);
    }

    fp_decimal_calculator::~fp_decimal_calculator() {
        mpz_clear(value);
    }

    fp_decimal_calculator::fp_decimal_calculator(const fp_decimal_calculator &other) : exp2(other.exp2) {
        mpz_init_set(value, other.value);
    }

    fp_decimal_calculator &fp_decimal_calculator::operator=(const fp_decimal_calculator &other) {
        if (this != &other) {
            mpz_set(value, other.value);
            exp2 = other.exp2;
        }
        return *this;
    }

    fp_decimal_calculator::fp_decimal_calculator(fp_decimal_calculator &&other) noexcept : exp2(other.exp2) {
        mpz_init(value);
        mpz_swap(value, other.value);
    }

    fp_decimal_calculator &fp_decimal_calculator::operator=(fp_decimal_calculator &&other) noexcept {
        if (this != &other) {
            mpz_swap(value, other.value);
            exp2 = other.exp2;
        }
        return *this;
    }

    void fp_decimal_calculator::fp_add(fp_decimal_calculator &out, const fp_decimal_calculator &a,
                                        const fp_decimal_calculator &b) noexcept {
        mpz_add(out.value, a.value, b.value);
    }

    void fp_decimal_calculator::fp_sub(fp_decimal_calculator &out, const fp_decimal_calculator &a,
                                        const fp_decimal_calculator &b) noexcept {
        mpz_sub(out.value, a.value, b.value);
    }

    void fp_decimal_calculator::fp_mul(fp_decimal_calculator &out, const fp_decimal_calculator &a,
                                        const fp_decimal_calculator &b) noexcept {
        mpz_t& temp = get_temp();
        mpz_mul(temp, a.value, b.value);
        mpz_tdiv_q_2exp(out.value, temp, static_cast<mp_bitcnt_t>(-a.exp2));
    }

    void fp_decimal_calculator::fp_div(fp_decimal_calculator &out, const fp_decimal_calculator &a,
                                        const fp_decimal_calculator &b) noexcept {
        mpz_t& temp = get_temp();
        const auto vbl = static_cast<int>(mpz_sizeinbase(b.value, 2));
        const int e = b.exp2 + vbl;

        mpz_mul_2exp(temp, a.value, static_cast<mp_bitcnt_t>(vbl));
        mpz_tdiv_q(temp, temp, b.value);

        if (e < 0) {
            mpz_mul_2exp(out.value, temp, static_cast<mp_bitcnt_t>(-e));
        } else {
            mpz_tdiv_q_2exp(out.value, temp, static_cast<mp_bitcnt_t>(e));
        }
    }

    void fp_decimal_calculator::fp_dbl(fp_decimal_calculator &out, const fp_decimal_calculator &target) noexcept {
        mpz_mul_2exp(out.value, target.value, 1);
    }

    void fp_decimal_calculator::fp_hlv(fp_decimal_calculator &out, const fp_decimal_calculator &target) noexcept {
        mpz_tdiv_q_2exp(out.value, target.value, 1);
    }

    void fp_decimal_calculator::fp_swap(fp_decimal_calculator &a, fp_decimal_calculator &b) noexcept {
        mpz_swap(a.value, b.value);
        std::swap(a.exp2, b.exp2);
    }

    void fp_decimal_calculator::negate(fp_decimal_calculator &target) noexcept {
        mpz_neg(target.value, target.value);
    }

    int fp_decimal_calculator::exp2ToExp10(const int e2) noexcept {
        return static_cast<int>(static_cast<double>(e2) * LOG10_2);
    }

    int fp_decimal_calculator::exp10ToExp2(const int precision) noexcept {
        return static_cast<int>(static_cast<double>(precision) / LOG10_2);
    }

    double fp_decimal_calculator::double_value() const noexcept {
        const int sgn = mpz_sgn(value);
        if (sgn == 0) return 0.0;

        mpz_t& temp = get_temp();
        mpz_abs(temp, value);

        const auto len = static_cast<int>(mpz_sizeinbase(temp, 2));
        const int shift = len - 53;

        if (shift < 0) {
            mpz_mul_2exp(temp, temp, static_cast<mp_bitcnt_t>(-shift));
        } else if (shift > 0) {
            mpz_tdiv_q_2exp(temp, temp, static_cast<mp_bitcnt_t>(shift));
        }

        uint64_t mantissa = 0;
        size_t cnt = 0;
        mpz_export(&mantissa, &cnt, -1, sizeof(mantissa), 0, 0, temp);

        const int fExp2 = exp2 + shift + 52;

        // Handle overflow
        if (fExp2 > 0x03ff) {
            return sgn > 0 ? std::numeric_limits<double>::infinity()
                           : -std::numeric_limits<double>::infinity();
        }

        // Handle denormals
        const int mantissa_shift = (fExp2 <= -0x03ff) ? (-0x03ff - fExp2 + 1) : 0;
        mantissa = (mantissa >> mantissa_shift) & 0x000fffffffffffffULL;

        const uint64_t exponent = (fExp2 <= -0x03ff)
            ? 0ULL
            : (0x3ff0000000000000ULL + (static_cast<uint64_t>(fExp2) << 52));

        const uint64_t sig = (sgn > 0) ? 0ULL : 0x8000000000000000ULL;

        return bit_cast_compat<double>(sig | exponent | mantissa);
    }

    void fp_decimal_calculator::double_exp_value(dex *result) const {
        const int sgn = mpz_sgn(value);
        if (sgn == 0) {
            dex::cpy(result, dex::ZERO);
            return;
        }

        mpz_t& temp = get_temp();
        mpz_abs(temp, value);

        const auto len = static_cast<int>(mpz_sizeinbase(temp, 2));
        const int shift = len - 53;

        if (shift < 0) {
            mpz_mul_2exp(temp, temp, static_cast<mp_bitcnt_t>(-shift));
        } else if (shift > 0) {
            mpz_tdiv_q_2exp(temp, temp, static_cast<mp_bitcnt_t>(shift));
        }

        uint64_t mantissa_bit = 0;
        size_t cnt = 0;
        mpz_export(&mantissa_bit, &cnt, -1, sizeof(mantissa_bit), 0, 0, temp);

        const int fExp2 = exp2 + shift + 52;
        const double mantissa = bit_cast_compat<double>(0x3ff0000000000000ULL | (mantissa_bit & 0x000fffffffffffffULL));

        dex::cpy(result, mantissa);
        dex::mul_2exp(result, *result, fExp2);
        if (sgn < 0) dex::neg(result);
    }

    void fp_decimal_calculator::setExp10(const int exp10) {
        const int newExp2 = exp10ToExp2(exp10);
        const int d_exp2 = exp2 - newExp2;

        if (d_exp2 < 0) {
            mpz_tdiv_q_2exp(value, value, static_cast<mp_bitcnt_t>(-d_exp2));
        } else if (d_exp2 > 0) {
            mpz_mul_2exp(value, value, static_cast<mp_bitcnt_t>(d_exp2));
        }
        exp2 = newExp2;
    }
}
// Vector expf and erff for the GNN, with the same results as glibc's scalar
// functions on a CPU with AVX2 and FMA, where glibc uses its FMA build of expf
// and erff calls that expf. Every lane does the glibc operations in the same
// order; checked against glibc 2.39 for all 2^32 inputs.
// expf: sysdeps/ieee754/flt-32/e_expf.c, erff: sysdeps/ieee754/flt-32/s_erff.c
//
// At first use the vector code is compared with the C library on probe values,
// including some that glibc does not round correctly. On any difference (another
// C library, or no AVX2 and FMA) the scalar functions are used instead.
#ifndef GNN_VEC_MATH_H
#define GNN_VEC_MATH_H

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

#if defined(GNN_VEC_MATH_IMPL) && defined(__x86_64__) && defined(__GNUC__) && !defined(__clang__)
#define GNN_VEC_MATH 1
#include <immintrin.h>
#endif

namespace gnn {
namespace vm {

void exp_n(float* p, size_t n);    // p[i] = std::exp(p[i])
// p[i] = 0.5f * p[i] * (1.0f + erff(p[i] * 0.70710678118654752f)): exact GELU, as
// torch.nn.functional.gelu (the tanh approximation differs by ~1e-3)
void gelu_n(float* p, size_t n);

#ifdef GNN_VEC_MATH_IMPL

#ifdef GNN_VEC_MATH

// glibc's __exp2f_data: 2^(i/32) table, shift, ln2 scaling and polynomial
static const uint64_t kExp2fTab[32] = {
0x3ff0000000000000, 0x3fefd9b0d3158574, 0x3fefb5586cf9890f, 0x3fef9301d0125b51,
0x3fef72b83c7d517b, 0x3fef54873168b9aa, 0x3fef387a6e756238, 0x3fef1e9df51fdee1,
0x3fef06fe0a31b715, 0x3feef1a7373aa9cb, 0x3feedea64c123422, 0x3feece086061892d,
0x3feebfdad5362a27, 0x3feeb42b569d4f82, 0x3feeab07dd485429, 0x3feea47eb03a5585,
0x3feea09e667f3bcd, 0x3fee9f75e8ec5f74, 0x3feea11473eb0187, 0x3feea589994cce13,
0x3feeace5422aa0db, 0x3feeb737b0cdc5e5, 0x3feec49182a3f090, 0x3feed503b23e255d,
0x3feee89f995ad3ad, 0x3feeff76f2fb5e47, 0x3fef199bdd85529c, 0x3fef3720dcef9069,
0x3fef5818dcfba487, 0x3fef7c97337b9b5f, 0x3fefa4afa2a490da, 0x3fefd0765b6e4540,
};
static const double kExpShift = 6755399441055744.0;
static const double kExpInvLn2N = 1.4426950408889634 * 32;
static const double kExpC0 = 0.055503615593415351 / 32 / 32 / 32;
static const double kExpC1 = 0.2402284522445722 / 32 / 32;
static const double kExpC2 = 0.69314718069162029 / 32;

// ---- AVX-512: 16 lanes, the expf core on two halves of 8 doubles ----
#pragma GCC push_options
#pragma GCC target("avx512f,avx2,fma")
namespace avx512 {

static inline __m256 exp_core8(__m256 x) {
    const __m512d xd = _mm512_cvtps_pd(x);
    __m512d kd = _mm512_fmadd_pd(_mm512_set1_pd(kExpInvLn2N), xd, _mm512_set1_pd(kExpShift));
    const __m512i ki = _mm512_castpd_si512(kd);
    kd = _mm512_sub_pd(kd, _mm512_set1_pd(kExpShift));
    const __m512d r = _mm512_fmsub_pd(_mm512_set1_pd(kExpInvLn2N), xd, kd);
    __m512i t = _mm512_i64gather_epi64(_mm512_and_si512(ki, _mm512_set1_epi64(31)),
                                       (const long long*)kExp2fTab, 8);
    t = _mm512_add_epi64(t, _mm512_slli_epi64(ki, 47));
    const __m512d z = _mm512_fmadd_pd(_mm512_set1_pd(kExpC0), r, _mm512_set1_pd(kExpC1));
    const __m512d r2 = _mm512_mul_pd(r, r);
    __m512d y = _mm512_fmadd_pd(r, _mm512_set1_pd(kExpC2), _mm512_set1_pd(1.0));
    y = _mm512_fmadd_pd(z, r2, y);
    return _mm512_cvtpd_ps(_mm512_mul_pd(y, _mm512_castsi512_pd(t)));
}
static inline __m512 exp_core(__m512 x) {
    const __m256 lo = exp_core8(_mm512_castps512_ps256(x));
    const __m256 hi = exp_core8(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(x), 1)));
    return _mm512_castpd_ps(_mm512_insertf64x4(_mm512_castps_pd(_mm512_castps256_ps512(lo)),
                                               _mm256_castps_pd(hi), 1));
}

#define W 16
#define F __m512
#define M __mmask16
#define SET1 _mm512_set1_ps
#define ADD _mm512_add_ps
#define SUB _mm512_sub_ps
#define MUL _mm512_mul_ps
#define DIV _mm512_div_ps
#define ZERO _mm512_setzero_ps()
#define ABSBITS(x) _mm512_and_si512(_mm512_castps_si512(x), _mm512_set1_epi32(0x7fffffff))
#define ABSBITS_TOP12(x) _mm512_and_si512(_mm512_srli_epi32(_mm512_castps_si512(x), 20), _mm512_set1_epi32(0x7ff))
#define ABS(x) _mm512_castsi512_ps(ABSBITS(x))
#define TRUNC12(x) _mm512_castsi512_ps(_mm512_and_si512(_mm512_castps_si512(x), _mm512_set1_epi32((int)0xfffff000)))
#define ILT(a, c) _mm512_cmplt_epi32_mask(a, _mm512_set1_epi32(c))
#define IGE(a, c) _mm512_cmpge_epi32_mask(a, _mm512_set1_epi32(c))
#define NEG(x) _mm512_cmplt_epi32_mask(_mm512_castps_si512(x), _mm512_setzero_si512())
#define MAND(a, b) (M)((a) & (b))
#define MOR(a, b) (M)((a) | (b))
#define MNOT(a) (M)(~(a))
#define ANY(m) ((m) != 0)
#define BITS(m) (unsigned)(m)
#define SEL(m, a, b) _mm512_mask_mov_ps(a, m, b)
#define LOAD _mm512_loadu_ps
#define STORE _mm512_storeu_ps
#define EXP_CORE exp_core
#include "GNNVecMath.inc"
#undef W
#undef F
#undef M
#undef SET1
#undef ADD
#undef SUB
#undef MUL
#undef DIV
#undef ZERO
#undef ABSBITS
#undef ABSBITS_TOP12
#undef ABS
#undef TRUNC12
#undef ILT
#undef IGE
#undef NEG
#undef MAND
#undef MOR
#undef MNOT
#undef ANY
#undef BITS
#undef SEL
#undef LOAD
#undef STORE
#undef EXP_CORE

} // namespace avx512
#pragma GCC pop_options

// ---- AVX2: 8 lanes, the expf core on two halves of 4 doubles ----
#pragma GCC push_options
#pragma GCC target("avx2,fma")
namespace avx2 {

static inline __m128 exp_core4(__m128 x) {
    const __m256d xd = _mm256_cvtps_pd(x);
    __m256d kd = _mm256_fmadd_pd(_mm256_set1_pd(kExpInvLn2N), xd, _mm256_set1_pd(kExpShift));
    const __m256i ki = _mm256_castpd_si256(kd);
    kd = _mm256_sub_pd(kd, _mm256_set1_pd(kExpShift));
    const __m256d r = _mm256_fmsub_pd(_mm256_set1_pd(kExpInvLn2N), xd, kd);
    __m256i t = _mm256_i64gather_epi64((const long long*)kExp2fTab,
                                       _mm256_and_si256(ki, _mm256_set1_epi64x(31)), 8);
    t = _mm256_add_epi64(t, _mm256_slli_epi64(ki, 47));
    const __m256d z = _mm256_fmadd_pd(_mm256_set1_pd(kExpC0), r, _mm256_set1_pd(kExpC1));
    const __m256d r2 = _mm256_mul_pd(r, r);
    __m256d y = _mm256_fmadd_pd(r, _mm256_set1_pd(kExpC2), _mm256_set1_pd(1.0));
    y = _mm256_fmadd_pd(z, r2, y);
    return _mm256_cvtpd_ps(_mm256_mul_pd(y, _mm256_castsi256_pd(t)));
}
static inline __m256 exp_core(__m256 x) {
    const __m128 lo = exp_core4(_mm256_castps256_ps128(x));
    const __m128 hi = exp_core4(_mm256_extractf128_ps(x, 1));
    return _mm256_insertf128_ps(_mm256_castps128_ps256(lo), hi, 1);
}

#define W 8
#define F __m256
#define M __m256i
#define SET1 _mm256_set1_ps
#define ADD _mm256_add_ps
#define SUB _mm256_sub_ps
#define MUL _mm256_mul_ps
#define DIV _mm256_div_ps
#define ZERO _mm256_setzero_ps()
#define ABSBITS(x) _mm256_and_si256(_mm256_castps_si256(x), _mm256_set1_epi32(0x7fffffff))
#define ABSBITS_TOP12(x) _mm256_and_si256(_mm256_srli_epi32(_mm256_castps_si256(x), 20), _mm256_set1_epi32(0x7ff))
#define ABS(x) _mm256_castsi256_ps(ABSBITS(x))
#define TRUNC12(x) _mm256_castsi256_ps(_mm256_and_si256(_mm256_castps_si256(x), _mm256_set1_epi32((int)0xfffff000)))
#define ILT(a, c) _mm256_cmpgt_epi32(_mm256_set1_epi32(c), a)
#define IGE(a, c) _mm256_xor_si256(ILT(a, c), _mm256_set1_epi32(-1))
#define NEG(x) _mm256_cmpgt_epi32(_mm256_setzero_si256(), _mm256_castps_si256(x))
#define MAND(a, b) _mm256_and_si256(a, b)
#define MOR(a, b) _mm256_or_si256(a, b)
#define MNOT(a) _mm256_xor_si256(a, _mm256_set1_epi32(-1))
#define ANY(m) (!_mm256_testz_si256(m, m))
#define BITS(m) (unsigned)_mm256_movemask_ps(_mm256_castsi256_ps(m))
#define SEL(m, a, b) _mm256_blendv_ps(a, b, _mm256_castsi256_ps(m))
#define LOAD _mm256_loadu_ps
#define STORE _mm256_storeu_ps
#define EXP_CORE exp_core
#include "GNNVecMath.inc"
#undef W
#undef F
#undef M
#undef SET1
#undef ADD
#undef SUB
#undef MUL
#undef DIV
#undef ZERO
#undef ABSBITS
#undef ABSBITS_TOP12
#undef ABS
#undef TRUNC12
#undef ILT
#undef IGE
#undef NEG
#undef MAND
#undef MOR
#undef MNOT
#undef ANY
#undef BITS
#undef SEL
#undef LOAD
#undef STORE
#undef EXP_CORE

} // namespace avx2
#pragma GCC pop_options

#endif // x86-64 GCC

namespace scalar {
void exp_n(float* p, size_t n) { for (size_t i = 0; i < n; ++i) p[i] = std::exp(p[i]); }
void erf_n(float* p, size_t n) { for (size_t i = 0; i < n; ++i) p[i] = std::erf(p[i]); }
void gelu_n(float* p, size_t n) {
    for (size_t i = 0; i < n; ++i) p[i] = 0.5f * p[i] * (1.0f + std::erf(p[i] * 0.70710678118654752f));
}
} // namespace scalar

struct Funcs { void (*exp_n)(float*, size_t); void (*erf_n)(float*, size_t); void (*gelu_n)(float*, size_t); };

// true if f gives the C library's results on the probe values
static bool same_as_libm(const Funcs& f) {
    // inputs where glibc's expf or erff is not correctly rounded, then a grid
    static const float probe[] = {
        0.0146949673f, -0.0744336545f, 0.0950903445f, -0.203090101f, 1.54697371f,
        -3.89493513f, 55.3368492f, 0.100504346f, -0.247672901f,
        7.4505806e-09f, -1.48021599e-07f, 2.07509771e-07f, -6.60304158e-06f, 7.55485235e-06f,
        -1.6834425e-05f, 2.44489111e-05f, -2.82561541e-05f, 0.000149188767f, -0.0002405626f,
        0.000533930899f, -0.0114057921f, 0.015304409f, -0.0493996553f, 0.0830827132f,
        -0.786486328f,
    };
    const size_t np = sizeof(probe) / sizeof(probe[0]);
    const size_t n = np + 2048;
    float in[n], a[n], b[n];
    std::memcpy(in, probe, sizeof(probe));
    for (size_t i = np; i < n; ++i) in[i] = -100.0f + 0.0977f * (float)(i - np) + 1.3e-4f * (float)(i % 97);
    std::memcpy(a, in, sizeof in); std::memcpy(b, in, sizeof in);
    f.exp_n(a, n); scalar::exp_n(b, n);
    if (std::memcmp(a, b, sizeof a)) return false;
    for (size_t i = np; i < n; ++i) in[i] *= 0.06f;   // erff over about -6 .. 6
    std::memcpy(a, in, sizeof in); std::memcpy(b, in, sizeof in);
    f.erf_n(a, n); scalar::erf_n(b, n);
    if (std::memcmp(a, b, sizeof a)) return false;
    std::memcpy(a, in, sizeof in); std::memcpy(b, in, sizeof in);
    f.gelu_n(a, n); scalar::gelu_n(b, n);
    return std::memcmp(a, b, sizeof a) == 0;
}

static Funcs choose() {
    const Funcs plain = {scalar::exp_n, scalar::erf_n, scalar::gelu_n};
#ifdef GNN_VEC_MATH
    __builtin_cpu_init();
    if (!__builtin_cpu_supports("avx2") || !__builtin_cpu_supports("fma")) return plain;
    const Funcs v512 = {avx512::exp_n, avx512::erf_n, avx512::gelu_n},
                v256 = {avx2::exp_n, avx2::erf_n, avx2::gelu_n};
    if (__builtin_cpu_supports("avx512f") && same_as_libm(v512)) return v512;
    if (same_as_libm(v256)) return v256;
#endif
    return plain;
}

static const Funcs& funcs() {
    static const Funcs f = choose();
    return f;
}

void exp_n(float* p, size_t n) { funcs().exp_n(p, n); }
void gelu_n(float* p, size_t n) { funcs().gelu_n(p, n); }

#endif // GNN_VEC_MATH_IMPL

} // namespace vm
} // namespace gnn

#endif

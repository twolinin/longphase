// Body of GNNVecMath.h, included once per instruction set with these macros:
//   W        lanes per vector
//   F, M     float vector and lane mask types
//   SET1, ADD, SUB, MUL, DIV, ZERO          float vector ops
//   ABSBITS(x) the bits of x without sign, as an integer vector
//   ILT(a, b), IGE(a, b)  signed compare of an integer vector with a constant: mask
//   MAND, MOR, MNOT, ANY, BITS(m)            masks
//   SEL(m, a, b)  b where m is set, else a
//   NEG(x)       lanes with the sign bit set: mask
//   LOAD, STORE  unaligned
//   EXP_CORE(x)  the expf main path on all lanes (FMA, in double)
// Each lane gets the operations of glibc's scalar code, in the same order.

static inline F expf_v(F x) {
    const M special = IGE(ABSBITS_TOP12(x), 0x42b);   // |x| >= 88 or NaN
    F out = EXP_CORE(x);
    if (ANY(special)) {
        alignas(64) float in[W], o[W];
        STORE(in, x); STORE(o, out);
        const unsigned b = BITS(special);
        for (int i = 0; i < W; ++i) if (b >> i & 1) o[i] = std::exp(in[i]);
        out = LOAD(o);
    }
    return out;
}

// c[0] + s*(c[1] + s*(... + s*c[n-1])), the order of the glibc source
template <int N>
static inline F horner_v(F s, const float (&c)[N]) {
    F acc = SET1(c[N - 1]);
    for (int i = N - 2; i >= 0; --i) acc = ADD(SET1(c[i]), MUL(s, acc));
    return acc;
}

static inline F erff_v(F x) {
    const F one = SET1(1.0f);
    const F ax = ABS(x);
    const M neg = NEG(x);
    // |x| < 2^-28, inf and NaN: scalar erff
    const M special = MOR(ILT(ABSBITS(x), 0x31800000), IGE(ABSBITS(x), 0x7f800000));
    const M mA = ILT(ABSBITS(x), 0x3f580000);                    // |x| < 0.84375
    const M mB = MAND(ILT(ABSBITS(x), 0x3fa00000), MNOT(mA));    // < 1.25
    const M mC = IGE(ABSBITS(x), 0x40c00000);                    // >= 6
    const M mD = MNOT(MOR(MOR(mA, mB), mC));
    F out = ZERO;
    if (ANY(mA)) {
        static const float pp[] = {1.2837916613e-01f, -3.2504209876e-01f, -2.8481749818e-02f,
                                   -5.7702702470e-03f, -2.3763017452e-05f};
        static const float qq[] = {1.0f, 3.9791721106e-01f, 6.5022252500e-02f, 5.0813062117e-03f,
                                   1.3249473704e-04f, -3.9602282413e-06f};
        const F z = MUL(x, x);
        const F y = DIV(horner_v(z, pp), horner_v(z, qq));
        out = SEL(mA, out, ADD(x, MUL(x, y)));
    }
    if (ANY(mB)) {
        static const float pa[] = {-2.3621185683e-03f, 4.1485610604e-01f, -3.7220788002e-01f,
                                   3.1834661961e-01f, -1.1089469492e-01f, 3.5478305072e-02f,
                                   -2.1663755178e-03f};
        static const float qa[] = {1.0f, 1.0642088205e-01f, 5.4039794207e-01f, 7.1828655899e-02f,
                                   1.2617121637e-01f, 1.3637083583e-02f, 1.1984500103e-02f};
        const F s = SUB(ax, one);
        const F erx = SET1(8.4506291151e-01f);
        const F pq = DIV(horner_v(s, pa), horner_v(s, qa));
        out = SEL(mB, out, SEL(neg, ADD(erx, pq), SUB(SUB(ZERO, erx), pq)));
    }
    if (ANY(mC)) {
        const F tiny = SET1(1e-30f);
        out = SEL(mC, out, SEL(neg, SUB(one, tiny), SUB(tiny, one)));
    }
    if (ANY(mD)) {
        static const float ra[] = {-9.8649440333e-03f, -6.9385856390e-01f, -1.0558626175e+01f,
                                   -6.2375331879e+01f, -1.6239666748e+02f, -1.8460508728e+02f,
                                   -8.1287437439e+01f, -9.8143291473e+00f};
        static const float sa[] = {1.0f, 1.9651271820e+01f, 1.3765776062e+02f, 4.3456588745e+02f,
                                   6.4538726807e+02f, 4.2900814819e+02f, 1.0863500214e+02f,
                                   6.5702495575e+00f, -6.0424413532e-02f};
        static const float rb[] = {-9.8649431020e-03f, -7.9928326607e-01f, -1.7757955551e+01f,
                                   -1.6063638306e+02f, -6.3756646729e+02f, -1.0250950928e+03f,
                                   -4.8351919556e+02f};
        static const float sb[] = {1.0f, 3.0338060379e+01f, 3.2579251099e+02f, 1.5367296143e+03f,
                                   3.1998581543e+03f, 2.5530502930e+03f, 4.7452853394e+02f,
                                   -2.2440952301e+01f};
        // the other lanes compute on 2.0, which stays on expf's main path
        const F xa = SEL(mD, SET1(2.0f), ax);
        const M ma = ILT(ABSBITS(xa), 0x4036DB6E);                 // |x| < 1/0.35
        const F s = DIV(one, MUL(xa, xa));
        const F R = SEL(ma, horner_v(s, rb), horner_v(s, ra));
        const F S = SEL(ma, horner_v(s, sb), horner_v(s, sa));
        const F z = TRUNC12(xa);                                   // low 12 bits cleared
        const F e1 = expf_v(SUB(MUL(SUB(ZERO, z), z), SET1(0.5625f)));
        const F e2 = expf_v(ADD(MUL(SUB(z, xa), ADD(z, xa)), DIV(R, S)));
        const F q = DIV(MUL(e1, e2), xa);
        out = SEL(mD, out, SEL(neg, SUB(one, q), SUB(q, one)));
    }
    if (ANY(special)) {
        alignas(64) float in[W], o[W];
        STORE(in, x); STORE(o, out);
        const unsigned b = BITS(special);
        for (int i = 0; i < W; ++i) if (b >> i & 1) o[i] = std::erf(in[i]);
        out = LOAD(o);
    }
    return out;
}

// 0.5f * x * (1.0f + erff(x * 0.70710678118654752f))
static inline F gelu_v(F v) {
    return MUL(MUL(SET1(0.5f), v), ADD(SET1(1.0f), erff_v(MUL(v, SET1(0.70710678118654752f)))));
}

// p[0..n) = f(p[0..n)), the last partial vector through a buffer padded with 1.0
template <F (*f)(F)>
static inline void apply_v(float* p, size_t n) {
    size_t i = 0;
    for (; i + W <= n; i += W) STORE(p + i, f(LOAD(p + i)));
    if (i < n) {
        alignas(64) float buf[W];
        for (int k = 0; k < W; ++k) buf[k] = 1.0f;
        std::memcpy(buf, p + i, (n - i) * sizeof(float));
        STORE(buf, f(LOAD(buf)));
        std::memcpy(p + i, buf, (n - i) * sizeof(float));
    }
}

void exp_n(float* p, size_t n) { apply_v<expf_v>(p, n); }
void erf_n(float* p, size_t n) { apply_v<erff_v>(p, n); }
void gelu_n(float* p, size_t n) { apply_v<gelu_v>(p, n); }

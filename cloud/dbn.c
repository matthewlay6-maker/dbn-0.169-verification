/* dbn.c -- certified (arb ball arithmetic) kernels for the Lambda <= 0.18 computation.
 *
 *   dbn block   Nm Np t0 y0                         -> rigorous lower bound for F_{Nm,Np}
 *   dbn barrier X0 y0 t0 ta tb level margin         -> rigorous barrier stepping on [ta, tb]
 *
 * Mathematics follows Polymath15 (arXiv:1904.12438) Section 8 exactly as implemented (float64)
 * in ../full_run.py (block_F) and ../barrier_gen.py; every quantity here is an arb/acb ball and every
 * inequality is decided on ball endpoints.  t0, y0 are passed as decimal strings and read exactly
 * (arb_set_str gives a ball containing the decimal).
 * Output: one JSON object per line on stdout.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <flint/arb.h>
#include <flint/acb.h>

/* mollifier prime set (env DBN_PRIMES, default 2,3,5) and constants (env DBN_TIGHT=1: tight constants
   sigma slack 1e-9 and |gamma| n^y <= e^{0.02 y}(1+1e-9) (n/N)^y, both from Polymath15 Thm 1.3 (20),(21);
   default: Polymath15's 0.001 and 1.03) */
static ulong G_P[8] = {2, 3, 5};
static int G_NP = 3;
static ulong DIVS[4096];
/* long mollifier from file (env DBN_MOLL: lines 'd lambda_d', lambda_1 = 1): any real coefficients are allowed */
static int G_LIST = 0;
static char *G_LAMSTR[4096];
static int NDIV = 0;
static ulong DPROD = 1;
static int G_TIGHT = 0;

static void setup_primes(void)
{
    const char *e = getenv("DBN_PRIMES");
    if (e && *e) {
        G_NP = 0;
        char buf[128]; strncpy(buf, e, 127); buf[127] = 0;
        for (char *tok = strtok(buf, ","); tok && G_NP < 8; tok = strtok(0, ",")) G_P[G_NP++] = strtoul(tok, 0, 10);
    }
    const char *tt = getenv("DBN_TIGHT");
    G_TIGHT = (tt && atoi(tt) == 1);
    const char *mf = getenv("DBN_MOLL");
    if (mf && *mf) {
        FILE *fh = fopen(mf, "r");
        if (!fh) { fprintf(stderr, "cannot open DBN_MOLL %s\n", mf); exit(3); }
        char line[512]; NDIV = 0; DPROD = 1;
        while (fgets(line, sizeof line, fh) && NDIV < 4096) {
            char ls[400]; unsigned long d;
            if (sscanf(line, "%lu %399s", &d, ls) != 2) continue;
            DIVS[NDIV] = d; G_LAMSTR[NDIV] = strdup(ls);
            if (d > DPROD) DPROD = d;
            NDIV++;
        }
        fclose(fh);
        if (NDIV == 0 || DIVS[0] != 1 || strcmp(G_LAMSTR[0], "1") != 0) { fprintf(stderr, "DBN_MOLL must start with '1 1'\n"); exit(3); }
        G_LIST = 1;
        return;
    }
    NDIV = 0; DPROD = 1;
    for (int i = 0; i < G_NP; i++) DPROD *= G_P[i];
    for (int mask = 0; mask < (1 << G_NP); mask++) {
        ulong d = 1;
        for (int i = 0; i < G_NP; i++) if (mask & (1 << i)) d *= G_P[i];
        DIVS[NDIV++] = d;
    }
}

static double ub(const arb_t x) { arf_t u; arf_init(u); arb_get_ubound_arf(u, x, 64); double d = arf_get_d(u, ARF_RND_UP); arf_clear(u); return d; }
static double lb(const arb_t x) { arf_t u; arf_init(u); arb_get_lbound_arf(u, x, 64); double d = arf_get_d(u, ARF_RND_DOWN); arf_clear(u); return d; }

/* ------------------------------------------------------------------ block */
/* block sum over n in [slo, shi): with slo = 2, shi = 30 Np + 1 this is the whole block.
   Partial runs print S_upper (rigorous upper bound of the partial sum); the parts are recombined in exact
   rational arithmetic by assemble_parts.py: F >= 1 - ca1_upper - sum S_upper - Z_upper. */
static int run_block(ulong Nm, ulong Np, const char *ts, const char *ys, ulong slo, ulong shi)
{
    slong prec = 80;
    arb_t t, y, lam[4096], alam[4096], c, sig, a1, r, tmp, tmp2, lm, b, a, w, u, v, total, one, ln;
    arb_init(t); arb_init(y); arb_init(c); arb_init(sig); arb_init(a1); arb_init(r);
    arb_init(tmp); arb_init(tmp2); arb_init(lm); arb_init(b); arb_init(a); arb_init(w);
    arb_init(u); arb_init(v); arb_init(total); arb_init(one); arb_init(ln);
    arb_set_str(t, ts, prec); arb_set_str(y, ys, prec); arb_one(one);
    for (int i = 0; i < NDIV; i++) {        /* lambda_d = prod_{p | d} (-b_p^t) */
        arb_init(lam[i]); arb_init(alam[i]); arb_one(lam[i]);
        if (G_LIST) { arb_set_str(lam[i], G_LAMSTR[i], prec); arb_abs(alam[i], lam[i]); continue; }
        for (int j = 0; j < G_NP; j++) if (DIVS[i] % G_P[j] == 0) {
            arb_log_ui(tmp, G_P[j], prec); arb_sqr(tmp, tmp, prec); arb_mul(tmp, tmp, t, prec);
            arb_mul_2exp_si(tmp, tmp, -2); arb_exp(tmp, tmp, prec); arb_neg(tmp, tmp);
            arb_mul(lam[i], lam[i], tmp, prec);
        }
        arb_abs(alam[i], lam[i]);
    }
    arb_log_ui(ln, Nm, prec);
    arb_mul(c, y, ln, prec); arb_neg(c, c); arb_exp(c, c, prec);             /* c = Nm^{-y0} */
    if (G_TIGHT) {                                                            /* alpha_1 = gcon */
        arb_set_str(tmp, "0.02", prec); arb_mul(tmp, tmp, y, prec); arb_exp(a1, tmp, prec);
        arb_set_str(tmp, "1.000000001", prec); arb_mul(a1, a1, tmp, prec);
    } else { arb_set_ui(a1, 103); arb_div_ui(a1, a1, 100, prec); }
    arb_mul(tmp, t, ln, prec); arb_mul_2exp_si(tmp, tmp, -1);                /* sigma = (1+y)/2 - 0.001 + (t/2) ln Nm */
    arb_add(sig, one, y, prec); arb_mul_2exp_si(sig, sig, -1); arb_add(sig, sig, tmp, prec);
    if (G_TIGHT) arb_set_str(tmp, "1e-9", prec); else { arb_set_ui(tmp, 1); arb_div_ui(tmp, tmp, 1000, prec); }
    arb_sub(sig, sig, tmp, prec);
    arb_mul(tmp, c, a1, prec); arb_sub(r, one, tmp, prec); arb_add(tmp2, one, tmp, prec); arb_div(r, r, tmp2, prec);
    ulong L = DPROD * Np, CH = 1UL << 20;
    if (slo < 2) slo = 2;
    if (shi == 0 || shi > L + 1) shi = L + 1;
    arb_ptr bc = _arb_vec_init(CH), ac = _arb_vec_init(CH), var = _arb_vec_init(CH);
    arb_zero(total);
    for (ulong s = slo; s < shi; s += CH) {
        ulong e = s + CH; if (e > shi) e = shi;
        ulong len = e - s;
        _arb_vec_zero(bc, len); _arb_vec_zero(ac, len); _arb_vec_zero(var, len);
        for (int i = 0; i < NDIV; i++) {
            ulong d = DIVS[i], mlo = (s + d - 1) / d, mhi = (e - 1) / d;
            if (mhi > Np) mhi = Np;
            for (ulong m = mlo; m <= mhi; m++) {
                arb_log_ui(lm, m, prec);
                arb_sqr(tmp, lm, prec); arb_mul(tmp, tmp, t, prec); arb_mul_2exp_si(tmp, tmp, -2);
                arb_exp(b, tmp, prec);                                       /* b_m = exp(t/4 ln^2 m) */
                arb_mul(tmp, y, lm, prec); arb_exp(tmp, tmp, prec);
                arb_mul(a, tmp, b, prec); arb_mul(a, a, a1, prec);           /* a_m = 1.03 m^y b_m */
                ulong k = m * d - s;
                if (m <= Nm) {
                    arb_addmul(bc + k, lam[i], b, prec);
                    arb_addmul(ac + k, lam[i], a, prec);
                } else {                                                     /* N-dependent: worst case */
                    arb_mul(tmp, c, a, prec); arb_add(tmp, tmp, b, prec);
                    arb_addmul(var + k, alam[i], tmp, prec);
                }
            }
        }
        for (ulong k = 0; k < len; k++) {
            arb_mul(tmp, c, ac + k, prec);
            arb_sub(u, bc + k, tmp, prec); arb_abs(u, u);
            arb_add(v, bc + k, tmp, prec); arb_abs(v, v); arb_mul(v, v, r, prec);
            arb_max(u, u, v, prec); arb_add(u, u, var + k, prec);
            if (arb_is_zero(u)) continue;
            arb_log_ui(lm, s + k, prec); arb_mul(lm, lm, sig, prec); arb_neg(lm, lm); arb_exp(w, lm, prec);
            arb_addmul(total, u, w, prec);
        }
    }
    /* Z (kappa term): Z <= 1.7 (N^{|kappa|}-1) sum_{n<=Np} n^{-q},  q = (1+y)/2 - 0.001,
       |kappa| <= t/(2(x_N - 6)),  x_N >= 4 pi Nm^2 - pi t/4 */
    arb_t kap, Z, q, pi;
    arb_init(kap); arb_init(Z); arb_init(q); arb_init(pi);
    arb_const_pi(pi, prec);
    arb_set_ui(tmp, Nm); arb_sqr(tmp, tmp, prec); arb_mul(tmp, tmp, pi, prec); arb_mul_2exp_si(tmp, tmp, 2);
    arb_mul(tmp2, pi, t, prec); arb_mul_2exp_si(tmp2, tmp2, -2); arb_sub(tmp, tmp, tmp2, prec);
    arb_sub_ui(tmp, tmp, 6, prec); arb_mul_2exp_si(tmp, tmp, 1); arb_div(kap, t, tmp, prec);
    arb_log_ui(lm, Np, prec); arb_mul(tmp, kap, lm, prec); arb_expm1(tmp, tmp, prec);
    arb_add(q, one, y, prec); arb_mul_2exp_si(q, q, -1); arb_set_ui(tmp2, 1); arb_div_ui(tmp2, tmp2, 1000, prec); arb_sub(q, q, tmp2, prec);
    arb_sub(tmp2, one, q, prec); arb_mul(Z, tmp2, lm, prec); arb_exp(Z, Z, prec); arb_div(Z, Z, tmp2, prec);  /* Np^{1-q}/(1-q) */
    arb_add_ui(Z, Z, 1, prec); arb_mul(Z, Z, tmp, prec);
    /* factor >= max|E| * gcon: 1.7 for Euler-5 (Polymath); 2 * a1 in general (max|E| < 2 for sigma >= 1.6, p <= 11) */
    if (G_LIST) {                       /* factor = sum |lambda_d| d^{-sigma} * gcon * 1.0001 >= max|M| * gcon */
        arb_t mh, dd; arb_init(mh); arb_init(dd); arb_zero(mh);
        for (int i = 0; i < NDIV; i++) { arb_log_ui(dd, DIVS[i], prec); arb_mul(dd, dd, sig, prec); arb_neg(dd, dd); arb_exp(dd, dd, prec); arb_addmul(mh, alam[i], dd, prec); }
        arb_mul(mh, mh, a1, prec); arb_set_str(dd, "1.0001", prec); arb_mul(mh, mh, dd, prec); arb_mul(Z, Z, mh, prec);
        arb_clear(mh); arb_clear(dd);
    }
    else if (G_TIGHT || G_NP != 3) { arb_mul_ui(tmp, a1, 2, prec); arb_mul(Z, Z, tmp, prec); }
    else { arb_set_ui(tmp, 17); arb_div_ui(tmp, tmp, 10, prec); arb_mul(Z, Z, tmp, prec); }
    /* F = 1 - c a1 - total - Z */
    arb_t F; arb_init(F);
    arb_mul(tmp, c, a1, prec); arb_sub(F, one, tmp, prec); arb_sub(F, F, total, prec); arb_sub(F, F, Z, prec);
    arb_mul(tmp, c, a1, prec);
    printf("{\"mode\":\"block\",\"Nminus\":%lu,\"Nplus\":%lu,\"slo\":%lu,\"shi\":%lu,\"full\":%s,\"F_lower\":%.12e,\"F_upper\":%.12e,"
           "\"D\":%lu,\"tight\":%d,\"S_upper\":%.17e,\"ca1_upper\":%.17e,\"Z_upper\":%.17e,\"sigma\":%.12f}\n",
           Nm, Np, slo, shi, (slo == 2 && shi == L + 1) ? "true" : "false", lb(F), ub(F), DPROD, G_TIGHT, ub(total), ub(tmp), ub(Z),
           arf_get_d(arb_midref(sig), ARF_RND_NEAR));
    fflush(stdout);
    return 0;
}

/* ---------------------------------------------------------------- barrier */
static slong P = 160;
static arb_t PI;

static void alpha(acb_t res, const acb_t s)      /* 1/(2s) + 1/(s-1) + (1/2) Log(s/(2 pi)) */
{
    acb_t a, b; acb_init(a); acb_init(b);
    acb_mul_2exp_si(a, s, 1); acb_inv(a, a, P);
    acb_sub_ui(b, s, 1, P); acb_inv(b, b, P); acb_add(a, a, b, P);
    acb_div_arb(b, s, PI, P); acb_mul_2exp_si(b, b, -1); acb_log(b, b, P); acb_mul_2exp_si(b, b, -1);
    acb_add(res, a, b, P);
    acb_clear(a); acb_clear(b);
}

static void logM0(acb_t res, const acb_t s)      /* Polymath15 eq. (7) */
{
    acb_t a, b; arb_t c; acb_init(a); acb_init(b); arb_init(c);
    acb_log(a, s, P);
    acb_sub_ui(b, s, 1, P); acb_log(b, b, P); acb_add(a, a, b, P);
    arb_log(c, PI, P); acb_mul_arb(b, s, c, P); acb_mul_2exp_si(b, b, -1); acb_sub(a, a, b, P);
    arb_mul_2exp_si(c, PI, 1); arb_sqrt(c, c, P); arb_mul_2exp_si(c, c, -4); arb_log(c, c, P); acb_add_arb(a, a, c, P);
    acb_mul_2exp_si(b, s, -1); acb_log(b, b, P);
    { acb_t h; arb_t hf; acb_init(h); arb_init(hf);            /* (s/2 - 1/2) Log(s/2) */
      acb_mul_2exp_si(h, s, -1); arb_set_d(hf, 0.5); acb_sub_arb(h, h, hf, P);
      acb_mul(b, b, h, P); acb_clear(h); arb_clear(hf); }
    acb_add(a, a, b, P);
    acb_mul_2exp_si(b, s, -1); acb_sub(res, a, b, P);
    acb_clear(a); acb_clear(b); arb_clear(c);
}

/* e_A + e_B + e_{C,0} upper bound at (x, y, t)  (Polymath15 Thm 1.3, Prop. "estimates" (iv)-(vi));
   decreasing in x and y (Polymath15 Lemma elem-lem (vi)).  drop_t != 0 discards the -t/16 log^2 term,
   so the bound then holds for all 0 <= t' <= t. */
static void e_bound(arb_t res, const arb_t x, const arb_t y, const arb_t t, int drop_t)
{
    arb_t lx, N, d1, s, eab, ec0, tmp, tmp2; acb_t z;
    arb_init(lx); arb_init(N); arb_init(d1); arb_init(s); arb_init(eab); arb_init(ec0); arb_init(tmp); arb_init(tmp2); acb_init(z);
    arb_div(lx, x, PI, P); arb_mul_2exp_si(lx, lx, -2);                       /* x/4pi */
    arb_div_ui(tmp, t, 16, P); arb_add(N, lx, tmp, P); arb_sqrt(N, N, P); arb_floor(N, N, P);
    arb_log(lx, lx, P);
    /* e_A + e_B <= (e^{d1} - 1)(1 + 1.021) sum_{n<=N} n^{-(1+y)/2} <= (e^{d1}-1)(1 + e^{0.02}) (N+1)^s / s, s = (1-y)/2 */
    arb_sqr(d1, lx, P); arb_mul(d1, d1, t, P); arb_mul(d1, d1, t, P); arb_div_ui(d1, d1, 16, P);
    arb_set_str(tmp, "0.626", P); arb_add(d1, d1, tmp, P);
    arb_set_str(tmp, "6.66", P); arb_sub(tmp, x, tmp, P); arb_div(d1, d1, tmp, P);
    arb_sub_ui(s, y, 1, P); arb_neg(s, s); arb_mul_2exp_si(s, s, -1);
    arb_expm1(eab, d1, P);
    arb_set_str(tmp, "0.02", P); arb_exp(tmp, tmp, P); arb_set_str(tmp2, "1.0000001", P); arb_mul(tmp, tmp, tmp2, P); arb_add_ui(tmp, tmp, 1, P);
    arb_mul(eab, eab, tmp, P);
    arb_add_ui(tmp, N, 1, P); arb_log(tmp, tmp, P); arb_mul(tmp, tmp, s, P); arb_exp(tmp, tmp, P); arb_div(tmp, tmp, s, P);
    arb_mul(eab, eab, tmp, P);
    /* e_{C,0} <= (x/4pi)^{-(1+y)/4} exp(-t/16 lx^2 + (3|lx + i pi/2| + 3.58)/(x - 8.52)) (1 + 1.24(3^y+3^-y)/(N-0.125) + 6.92/(x-6.66)) */
    arb_add_ui(tmp, y, 1, P); arb_mul(tmp, tmp, lx, P); arb_mul_2exp_si(tmp, tmp, -2); arb_neg(ec0, tmp);
    if (!drop_t) { arb_sqr(tmp, lx, P); arb_mul(tmp, tmp, t, P); arb_div_ui(tmp, tmp, 16, P); arb_sub(ec0, ec0, tmp, P); }
    arb_mul_2exp_si(tmp2, PI, -1); acb_set_arb_arb(z, lx, tmp2); acb_abs(tmp, z, P); arb_mul_ui(tmp, tmp, 3, P);
    arb_set_str(tmp2, "3.58", P); arb_add(tmp, tmp, tmp2, P);
    arb_set_str(tmp2, "8.52", P); arb_sub(tmp2, x, tmp2, P); arb_div(tmp, tmp, tmp2, P); arb_add(ec0, ec0, tmp, P);
    arb_exp(ec0, ec0, P);
    arb_set_ui(tmp, 3); arb_pow(tmp, tmp, y, P); arb_inv(tmp2, tmp, P); arb_add(tmp, tmp, tmp2, P);
    arb_set_str(tmp2, "1.24", P); arb_mul(tmp, tmp, tmp2, P);
    arb_set_str(tmp2, "0.125", P); arb_sub(tmp2, N, tmp2, P); arb_div(tmp, tmp, tmp2, P); arb_add_ui(tmp, tmp, 1, P);
    arb_set_str(tmp2, "6.66", P); arb_sub(tmp2, x, tmp2, P); arb_set_str(s, "6.92", P); arb_div(tmp2, s, tmp2, P); arb_add(tmp, tmp, tmp2, P);
    arb_mul(ec0, ec0, tmp, P);
    arb_add(res, eab, ec0, P);
    arb_clear(lx); arb_clear(N); arb_clear(d1); arb_clear(s); arb_clear(eab); arb_clear(ec0); arb_clear(tmp); arb_clear(tmp2); acb_clear(z);
}

/* mollifier E_{t,5} on Re s >= sig: hi >= |E|, lo <= |E|, ar >= |arg E| */
static void mollifier(arb_t hi, arb_t lo, arb_t ar, const arb_t t, const arb_t sig)
{
    arb_t v, tmp; arb_init(v); arb_init(tmp);
    arb_one(hi); arb_one(lo); arb_zero(ar);
    if (G_LIST) {   /* M = 1 + R, |R| <= rho = sum_{d>=2} |lambda_d| d^{-sig}: |M| in [1-rho, 1+rho], |arg M| <= asin(rho) */
        arb_t rho, l; arb_init(rho); arb_init(l); arb_zero(rho);
        for (int i = 1; i < NDIV; i++) {
            arb_set_str(l, G_LAMSTR[i], P); arb_abs(l, l);
            arb_log_ui(tmp, DIVS[i], P); arb_mul(tmp, tmp, sig, P); arb_neg(tmp, tmp); arb_exp(tmp, tmp, P);
            arb_addmul(rho, l, tmp, P);
        }
        arb_add_ui(hi, rho, 1, P); arb_sub_ui(lo, rho, 1, P); arb_neg(lo, lo);
        arb_one(tmp);
        if (arb_lt(rho, tmp)) arb_asin(ar, rho, P); else arb_set(ar, PI);
        arb_clear(rho); arb_clear(l); arb_clear(v); arb_clear(tmp);
        return;
    }
    for (int j = 0; j < G_NP; j++) {
        arb_log_ui(tmp, G_P[j], P); arb_sqr(v, tmp, P); arb_mul(v, v, t, P); arb_mul_2exp_si(v, v, -2);
        arb_submul(v, sig, tmp, P); arb_exp(v, v, P);                          /* b_p p^{-sig} */
        arb_add_ui(tmp, v, 1, P); arb_mul(hi, hi, tmp, P);
        arb_sub_ui(tmp, v, 1, P); arb_neg(tmp, tmp); arb_mul(lo, lo, tmp, P);
        arb_asin(tmp, v, P); arb_add(ar, ar, tmp, P);
    }
    arb_clear(v); arb_clear(tmp);
}

/* lower bound for dist(z, (-inf, 0]) over all z with |z| >= rr and |arg z| <= ph; -1 if ph may reach pi */
static double dist_lb(const arb_t rr, const arb_t ph)
{
    arb_t d, tmp; arb_init(d); arb_init(tmp);
    arb_sub(tmp, PI, ph, P);
    if (!arb_is_positive(tmp)) { arb_clear(d); arb_clear(tmp); return -1.0; }
    arb_mul_2exp_si(tmp, PI, -1);
    if (arb_le(ph, tmp)) arb_set(d, rr);
    else { arb_sub(tmp, PI, ph, P); arb_sin(tmp, tmp, P); arb_mul(d, rr, tmp, P); }
    double out = lb(d);
    arb_clear(d); arb_clear(tmp);
    return out;
}

/* sig_min = (1+y)/2 + (t/2) log N0 - 1e-6 and need = E_hi * e(X, y, t) for region (ii) */
static void region_ii_budget(arb_t Ehi, arb_t Elo, arb_t Ear, arb_t need, const arb_t x, const arb_t y, const arb_t t, ulong N0)
{
    arb_t sig, tmp, e; arb_init(sig); arb_init(tmp); arb_init(e);
    arb_add_ui(sig, y, 1, P); arb_mul_2exp_si(sig, sig, -1);
    arb_log_ui(tmp, N0, P); arb_mul(tmp, tmp, t, P); arb_mul_2exp_si(tmp, tmp, -1); arb_add(sig, sig, tmp, P);
    arb_set_str(tmp, "1e-6", P); arb_sub(sig, sig, tmp, P);
    mollifier(Ehi, Elo, Ear, t, sig);
    e_bound(e, x, y, t, 0); arb_mul(need, Ehi, e, P);
    arb_clear(sig); arb_clear(tmp); arb_clear(e);
}

static int LEFT_ONLY = 0;

typedef struct { arb_t d, y; } pt_t;

static int run_barrier(const char *X0s, const char *ys, const char *t0s, double ta, double tb, int level, const char *margins)
{
    arb_init(PI); arb_const_pi(PI, P);
    arb_t X0, y0, YC, tmp, tmp2, x, lx, margin;
    arb_init(X0); arb_init(y0); arb_init(YC); arb_init(tmp); arb_init(tmp2); arb_init(x); arb_init(lx); arb_init(margin);
    arb_set_str(X0, X0s, P); arb_set_str(y0, ys, P); arb_set_str(margin, margins, P);
    arb_add_ui(YC, y0, 1, P); arb_mul_2exp_si(YC, YC, -1);
    /* N constant on the barrier: floor(sqrt((X0-1/2)/4pi)) == floor(sqrt((X0+1/2)/4pi + t0/16)) */
    arb_t t0; arb_init(t0); arb_set_str(t0, t0s, P);
    arb_set_d(tmp, 0.5); arb_sub(x, X0, tmp, P); arb_div(tmp, x, PI, P); arb_mul_2exp_si(tmp, tmp, -2); arb_sqrt(tmp, tmp, P);
    arb_set_d(tmp2, 0.5); arb_add(tmp2, X0, tmp2, P); arb_div(tmp2, tmp2, PI, P); arb_mul_2exp_si(tmp2, tmp2, -2);
    { arb_t q; arb_init(q); arb_div_ui(q, t0, 16, P); arb_add(tmp2, tmp2, q, P); arb_clear(q); }
    arb_sqrt(tmp2, tmp2, P);
    fmpz_t f1, f2; fmpz_init(f1); fmpz_init(f2);
    arb_floor(tmp, tmp, P); arb_floor(tmp2, tmp2, P);
    if (!arb_get_unique_fmpz(f1, tmp) || !arb_get_unique_fmpz(f2, tmp2) || !fmpz_equal(f1, f2)) {
        printf("{\"mode\":\"barrier\",\"error\":\"N not constant\"}\n"); return 1; }
    ulong N = fmpz_get_ui(f1);
    /* x_low = X0 - 1/2 used in derivative bounds; check N^2 <= x_low/4pi */
    arb_set_d(tmp, 0.5); arb_sub(x, X0, tmp, P);
    arb_div(lx, x, PI, P); arb_mul_2exp_si(lx, lx, -2);
    arb_set_ui(tmp, N); arb_sqr(tmp, tmp, P);
    if (!arb_le(tmp, lx)) { printf("{\"mode\":\"barrier\",\"error\":\"N^2 > x/4pi\"}\n"); return 1; }
    arb_log(lx, lx, P);                                   /* lx = log(x/4pi) */

    /* per-n data */
    arb_ptr ln = _arb_vec_init(N + 1);
    acb_ptr ph = _acb_vec_init(N + 1);
    arb_t half; arb_init(half);
    for (ulong n = 1; n <= N; n++) {
        arb_log_ui(ln + n, n, P);
        arb_mul(tmp, X0, ln + n, P); arb_mul_2exp_si(tmp, tmp, -1);   /* (X0/2) ln n */
        arb_sin_cos(acb_imagref(ph + n), acb_realref(ph + n), tmp, P);
    }
    /* mesh (counter-clockwise), spacing <= 2^-level */
    ulong M = 1UL << level;
    arb_t hgt; arb_init(hgt); arb_set_ui(hgt, 1); arb_sub(hgt, hgt, y0, P);
    ulong nb = M, ny = (ulong) ceil((1.0 - atof(ys)) * M);
    ulong npts = 2 * nb + 2 * ny;
    arb_ptr pd = _arb_vec_init(npts), py = _arb_vec_init(npts);
    ulong j = 0;
    for (ulong k = 0; k < nb; k++, j++) { arb_set_si(pd + j, 2 * (slong) k - (slong) nb); arb_div_ui(pd + j, pd + j, 2 * nb, P); arb_set(py + j, y0); }
    for (ulong k = 0; k < ny; k++, j++) { arb_set_d(pd + j, 0.5); arb_mul_ui(tmp, hgt, k, P); arb_div_ui(tmp, tmp, ny, P); arb_add(py + j, y0, tmp, P); }
    for (ulong k = 0; k < nb; k++, j++) { arb_set_si(pd + j, (slong) nb - 2 * (slong) k); arb_div_ui(pd + j, pd + j, 2 * nb, P); arb_one(py + j); }
    for (ulong k = 0; k < ny; k++, j++) { arb_set_d(pd + j, -0.5); arb_mul_ui(tmp, hgt, k, P); arb_div_ui(tmp, tmp, ny, P); arb_one(py + j); arb_sub(py + j, py + j, tmp, P); }
    double h0 = 1.0 / (double) M;
    /* gamma parts per point: G0 = M0(s)/M0(1-s), C = (alpha(s)^2 - alpha(1-s)^2)/4, s = (1-y+ix)/2; w = (i d - (y-YC))/2 */
    acb_ptr G0 = _acb_vec_init(npts), CC = _acb_vec_init(npts), W = _acb_vec_init(npts);
    acb_t s, s1, A, B; acb_init(s); acb_init(s1); acb_init(A); acb_init(B);
    for (j = 0; j < npts; j++) {
        arb_add(tmp, X0, pd + j, P);
        arb_sub_ui(tmp2, py + j, 1, P); arb_neg(tmp2, tmp2);
        acb_set_arb_arb(s, tmp2, tmp); acb_mul_2exp_si(s, s, -1);
        acb_one(s1); acb_sub(s1, s1, s, P);
        logM0(A, s); logM0(B, s1); acb_sub(A, A, B, P); acb_exp(G0 + j, A, P);
        alpha(A, s); alpha(B, s1); acb_sqr(A, A, P); acb_sqr(B, B, P); acb_sub(A, A, B, P); acb_mul_2exp_si(CC + j, A, -2);
        arb_sub(tmp, py + j, YC, P); arb_neg(tmp, tmp);
        acb_set_arb_arb(W + j, tmp, pd + j); acb_mul_2exp_si(W + j, W + j, -1);
    }
    fprintf(stderr, "barrier: N=%lu points=%lu precompute done\n", N, npts);

    const int K = 40;
    acb_ptr M1 = _acb_vec_init(K + 1), M2 = _acb_vec_init(K + 1);
    acb_ptr fv = _acb_vec_init(npts);
    arb_t Dz, Dt, sig, wn, big, q, g, abs1, abs2, tailrel, allow, L, fact, mn, spare, xm6, e02;
    arb_init(Dz); arb_init(Dt); arb_init(sig); arb_init(wn); arb_init(big); arb_init(q); arb_init(g);
    arb_init(abs1); arb_init(abs2); arb_init(tailrel); arb_init(allow); arb_init(L); arb_init(fact);
    arb_init(mn); arb_init(spare); arb_init(xm6); arb_init(e02);
    arb_sub_ui(xm6, x, 6, P);
    arb_set_d(e02, 0.02); arb_exp(e02, e02, P);
    /* big = log(sqrt(4 + x^2)/4pi) + pi + 3/x */
    arb_sqr(big, x, P); arb_add_ui(big, big, 4, P); arb_sqrt(big, big, P); arb_div(big, big, PI, P); arb_mul_2exp_si(big, big, -2); arb_log(big, big, P);
    arb_add(big, big, PI, P); arb_set_ui(tmp, 3); arb_div(tmp, tmp, x, P); arb_add(big, big, tmp, P);
    /* Taylor tail per unit |a_n| n^{0.34}: (0.34 L)^{K+1}/(K+1)! (|w| <= 0.34 checked below) */
    arb_log_ui(L, N, P); arb_set_d(tmp, 0.34); arb_mul(tmp, tmp, L, P); arb_pow_ui(tailrel, tmp, K + 1, P);
    arb_fac_ui(fact, K + 1, P); arb_div(tailrel, tailrel, fact, P);
    for (j = 0; j < npts; j++) { acb_abs(tmp, W + j, P); arb_set_d(tmp2, 0.34); if (!arb_le(tmp, tmp2)) { printf("{\"error\":\"|w|>0.34\"}\n"); return 1; } }

    double t = ta;
    long steps = 0; int ok = 1;
    double minf_all = 1e300;
    arb_t T, b, ex, re1, re2; arb_init(T); arb_init(b); arb_init(ex); arb_init(re1); arb_init(re2);
    acb_t a1, a2, xc, z, c1, c2, pw, fz, gm, acc1, acc2, wk;
    acb_init(a1); acb_init(a2); acb_init(xc); acb_init(z); acb_init(c1); acb_init(c2); acb_init(pw); acb_init(fz);
    acb_init(gm); acb_init(acc1); acb_init(acc2); acb_init(wk);
    arb_ptr pk = _arb_vec_init(K + 1);
    while (t < tb) {
        arb_set_d(T, t);
        /* derivative bounds (Polymath15 Lemma, barrier section), valid for all t' >= t */
        arb_mul(sig, T, lx, P); arb_mul_2exp_si(sig, sig, -2);
        arb_add_ui(tmp, y0, 1, P); arb_mul_2exp_si(tmp, tmp, -1); arb_add(sig, sig, tmp, P);
        arb_set_d(tmp, 1e-18); arb_sub(sig, sig, tmp, P);
        arb_mul_ui(g, e02, 1000001, P); arb_div_ui(g, g, 1000000, P);
        arb_zero(Dz); arb_zero(Dt);
        arb_t s1z, s2z, s1t, s2t, tq, cc; arb_init(s1z); arb_init(s2z); arb_init(s1t); arb_init(s2t); arb_init(tq); arb_init(cc);
        arb_mul_2exp_si(tq, xm6, 2); arb_div(tq, T, tq, P);                     /* t/(4(x-6)) */
        /* cc = 0.25 (pi/2 + 8/(x-6)) (lx + 8/(x-6)) */
        arb_set_ui(tmp, 8); arb_div(tmp, tmp, xm6, P);
        arb_mul_2exp_si(tmp2, PI, -1); arb_add(tmp2, tmp2, tmp, P); arb_add(cc, lx, tmp, P); arb_mul(cc, cc, tmp2, P); arb_mul_2exp_si(cc, cc, -2);
        for (ulong n = 1; n <= N; n++) {
            arb_sqr(tmp, ln + n, P); arb_mul(tmp, tmp, T, P); arb_mul_2exp_si(tmp, tmp, -2);
            arb_submul(tmp, sig, ln + n, P); arb_exp(wn, tmp, P);              /* b_n / n^{sigma} */
            /* Dz terms */
            arb_mul_2exp_si(tmp, ln + n, -1); arb_addmul(tmp, tq, ln + n, P); arb_addmul(s1z, wn, tmp, P);
            arb_mul(tmp, tq, ln + n, P); arb_set_d(tmp2, 0.5); arb_add(tmp2, tmp2, tq, P); arb_addmul(tmp, big, tmp2, P); arb_addmul(s2z, wn, tmp, P);
            /* Dt terms: q = 0.25 ln n ln(x/4pi n) + pi/8 ln n + 2 ln n/(x-6) */
            arb_sub(tmp, lx, ln + n, P); arb_mul(q, tmp, ln + n, P); arb_mul_2exp_si(q, q, -2);
            arb_mul(tmp, PI, ln + n, P); arb_mul_2exp_si(tmp, tmp, -3); arb_add(q, q, tmp, P);
            arb_mul_2exp_si(tmp, ln + n, 1); arb_div(tmp, tmp, xm6, P); arb_add(q, q, tmp, P);
            arb_addmul(s1t, wn, q, P);
            arb_add(tmp, q, cc, P); arb_addmul(s2t, wn, tmp, P);
        }
        arb_mul(tmp, g, s2z, P); arb_add(Dz, s1z, tmp, P);
        arb_mul(tmp, g, s2t, P); arb_add(Dt, s1t, tmp, P);
        double Dzu = ub(Dz), Dtu = ub(Dt);
        arb_clear(s1z); arb_clear(s2z); arb_clear(s1t); arb_clear(s2t); arb_clear(tq); arb_clear(cc);

        /* moments at the centre */
        arb_t yc1; arb_init(yc1);
        arb_add_ui(yc1, YC, 1, P);
        acb_set_arb_arb(xc, yc1, X0); acb_conj(xc, xc); acb_mul_2exp_si(xc, xc, -1);  /* (1+YC - i X0)/2 */
        alpha(a1, xc); acb_mul_arb(a1, a1, T, P); acb_mul_2exp_si(a1, a1, -1);
        arb_sub_ui(yc1, YC, 1, P); arb_neg(yc1, yc1);
        acb_set_arb_arb(xc, yc1, X0); acb_mul_2exp_si(xc, xc, -1);                    /* (1-YC + i X0)/2 */
        alpha(a2, xc); acb_mul_arb(a2, a2, T, P); acb_mul_2exp_si(a2, a2, -1);
        _acb_vec_zero(M1, K + 1); _acb_vec_zero(M2, K + 1); arb_zero(abs1); arb_zero(abs2);
        arb_add_ui(re1, YC, 1, P); arb_mul_2exp_si(re1, re1, -1);                     /* (1+YC)/2 */
        arb_sub_ui(re2, YC, 1, P); arb_neg(re2, re2); arb_mul_2exp_si(re2, re2, -1);  /* (1-YC)/2 */
        for (ulong n = 1; n <= N; n++) {
            arb_sqr(tmp, ln + n, P); arb_mul(tmp, tmp, T, P); arb_mul_2exp_si(tmp, tmp, -2); arb_exp(b, tmp, P);
            acb_add_arb(z, a1, re1, P); acb_mul_arb(z, z, ln + n, P); acb_neg(z, z); acb_exp(c1, z, P);   /* n^{-((1+YC)/2 + a1)} */
            acb_mul_arb(c1, c1, b, P); acb_mul(c1, c1, ph + n, P);
            acb_add_arb(z, a2, re2, P); acb_mul_arb(z, z, ln + n, P); acb_neg(z, z); acb_exp(c2, z, P);   /* n^{-((1-YC)/2 + a2)} */
            acb_mul_arb(c2, c2, b, P); acb_conj(pw, ph + n); acb_mul(c2, c2, pw, P);
            arb_one(pk + 0);
            for (int k = 1; k <= K; k++) arb_mul(pk + k, pk + k - 1, ln + n, P);
            for (int k = 0; k <= K; k++) { acb_addmul_arb(M1 + k, c1, pk + k, P); acb_addmul_arb(M2 + k, c2, pk + k, P); }
            arb_set_d(tmp2, 0.34); arb_mul(tmp2, tmp2, ln + n, P); arb_exp(tmp2, tmp2, P);
            acb_abs(tmp, c1, P); arb_addmul(abs1, tmp, tmp2, P);
            acb_abs(tmp, c2, P); arb_addmul(abs2, tmp, tmp2, P);
        }
        arb_clear(yc1);
        /* allowance: Taylor tail (rigorous) + alpha/kappa variation 2e-11 relative (see barrier_gen.py) */
        arb_mul(tmp, e02, abs2, P); arb_add(tmp, tmp, abs1, P);
        arb_set_d(tmp2, 2e-11); arb_add(tmp2, tmp2, tailrel, P); arb_mul(allow, tmp, tmp2, P);
        double allowu = ub(allow);

        if (LEFT_ONLY) {
            /* left edge of region (ii): x = X0 - 1/2, y0 <= y <= 1, t = t0 (the time passed as ta).
               Every boundary point is within h0/2 of a mesh point, so |f - f(z_j)| <= slack. */
            arb_t fmin, amax, slack, rr, ph, Ehi, Elo, Ear, need, ab, ar;
            arb_init(fmin); arb_init(amax); arb_init(slack); arb_init(rr); arb_init(ph);
            arb_init(Ehi); arb_init(Elo); arb_init(Ear); arb_init(need); arb_init(ab); arb_init(ar);
            double fminl = 1e300, amaxu = 0;
            ulong cnt = 0;
            for (j = 0; j < npts; j++) {
                if (!(j == 0 || j >= 2 * nb + ny)) continue;          /* left side incl. bottom-left corner */
                acb_zero(acc1); acb_zero(acc2); acb_neg(wk, W + j);
                for (int k = K; k >= 0; k--) {
                    acb_mul(acc1, acc1, W + j, P); acb_div_ui(acc1, acc1, k + 1, P); acb_add(acc1, acc1, M1 + k, P);
                    acb_mul(acc2, acc2, wk, P); acb_div_ui(acc2, acc2, k + 1, P); acb_add(acc2, acc2, M2 + k, P);
                }
                acb_mul_arb(gm, CC + j, T, P); acb_exp(gm, gm, P); acb_mul(gm, gm, G0 + j, P);
                acb_mul(fz, gm, acc2, P); acb_add(fz, fz, acc1, P);
                acb_abs(ab, fz, P); acb_arg(ar, fz, P); arb_abs(ar, ar);
                double l = lb(ab), a = ub(ar);
                if (l < fminl) fminl = l;
                if (a > amaxu) amaxu = a;
                cnt++;
            }
            /* slack = D_z h0/2 + allowance + D_t * 1e-15 (t0 rounded to double) */
            arb_set_d(slack, Dzu * h0 / 2 + allowu + Dtu * 1e-15);
            arb_set_d(fmin, fminl); arb_set_d(amax, amaxu);
            arb_sub(rr, fmin, slack, P);
            arb_div(ph, slack, fmin, P); arb_asin(ph, ph, P); arb_add(ph, ph, amax, P);
            region_ii_budget(Ehi, Elo, Ear, need, x, y0, T, N);
            arb_add(ph, ph, Ear, P); arb_mul(rr, rr, Elo, P);
            double d = dist_lb(rr, ph);
            int okL = d > ub(need);
            printf("{\"mode\":\"leftedge\",\"t\":%.17g,\"points\":%lu,\"min_abs_f\":%.10e,\"max_abs_arg\":%.10e,\"slack\":%.6e,"
                   "\"dist\":%.10e,\"need\":%.6e,\"ok\":%s}\n", t, cnt, fminl, amaxu, ub(slack), d, ub(need), okL ? "true" : "false");
            fflush(stdout);
            return okL ? 0 : 2;
        }

        /* evaluate on a coarse sub-mesh to pick the spacing */
        double mest = 1e300;
        ulong stride = 64;
        for (int pass = 0; pass < 2; pass++) {
            ulong step = (pass == 0) ? stride : 0;
            if (pass == 1) { step = stride; while (step > 1 && Dzu * (h0 * step) / 2 > 0.2 * mest) step /= 2; }
            double mloc = 1e300; double windsum = 0; int badarg = 0;
            ulong cnt = 0; acb_t prev; acb_init(prev);
            for (j = 0; j < npts; j += step) {
                /* f = sum_k M1[k] w^k/k! + gamma sum_k M2[k] (-w)^k/k!; Horner: acc = acc*w/(k+1) + M_k */
                acb_zero(acc1); acb_zero(acc2);
                acb_neg(wk, W + j);
                for (int k = K; k >= 0; k--) {
                    acb_mul(acc1, acc1, W + j, P); acb_div_ui(acc1, acc1, k + 1, P); acb_add(acc1, acc1, M1 + k, P);
                    acb_mul(acc2, acc2, wk, P); acb_div_ui(acc2, acc2, k + 1, P); acb_add(acc2, acc2, M2 + k, P);
                }
                acb_mul_arb(gm, CC + j, T, P); acb_exp(gm, gm, P); acb_mul(gm, gm, G0 + j, P);
                acb_mul(fz, gm, acc2, P); acb_add(fz, fz, acc1, P);
                acb_abs(tmp, fz, P); double lo = lb(tmp); if (lo < mloc) mloc = lo;
                if (pass == 1) {
                    acb_set(fv + cnt, fz);
                    if (cnt > 0) {
                        acb_div(z, fz, prev, P); acb_arg(tmp, z, P);
                        double a = ub(tmp), al = lb(tmp);
                        if (a > 1.5 || al < -1.5) badarg = 1;
                        windsum += arf_get_d(arb_midref(tmp), ARF_RND_NEAR);
                    }
                    acb_set(prev, fz); cnt++;
                }
            }
            if (pass == 1) {
                acb_div(z, fv + 0, prev, P); acb_arg(tmp, z, P);
                if (ub(tmp) > 1.5 || lb(tmp) < -1.5) badarg = 1;
                windsum += arf_get_d(arb_midref(tmp), ARF_RND_NEAR);
                long wind = lround(windsum / (2 * M_PI));
                double h = h0 * step;
                arb_set_d(mn, mloc);
                arb_set_d(tmp, Dzu * h / 2); arb_sub(spare, mn, tmp, P); arb_sub(spare, spare, margin, P);
                arb_set_d(tmp, allowu); arb_sub(spare, spare, tmp, P);
                double sp = lb(spare);
                if (mloc < minf_all) minf_all = mloc;
                if (sp <= 0 || wind != 0 || badarg) {
                    printf("{\"mode\":\"barrier\",\"fail\":true,\"t\":%.17g,\"min_abs_f\":%.6e,\"winding\":%ld,\"badarg\":%d}\n", t, mloc, wind, badarg);
                    ok = 0; break;
                }
                double dt = sp / Dtu * 0.999999;
                if (steps % 25 == 0)
                    printf("{\"mode\":\"barrier_step\",\"t\":%.17g,\"dt\":%.6e,\"points\":%lu,\"min_abs_f\":%.6e,\"Dz\":%.4e,\"Dt\":%.4e,\"allow\":%.3e}\n",
                           t, dt, cnt, mloc, Dzu, Dtu, allowu);
                fflush(stdout);
                t = t + dt; steps++;
            } else mest = mloc;
            acb_clear(prev);
        }
        if (!ok) break;
    }
    printf("{\"mode\":\"barrier\",\"ta\":%.17g,\"tb\":%.17g,\"t_reached\":%.17g,\"ok\":%s,\"steps\":%ld,\"min_abs_f\":%.6e,\"N\":%lu,\"points_fine\":%lu}\n",
           ta, tb, t, (ok && t >= tb) ? "true" : "false", steps, minf_all, N, npts);
    fflush(stdout);
    return ok ? 0 : 2;
}

/* left edge of region (ii) at t = t0 (full fine mesh), see LEFT_ONLY block in run_barrier */
static int run_leftedge(const char *X0s, const char *ys, const char *t0s, int level)
{
    LEFT_ONLY = 1;
    double t0 = atof(t0s);
    return run_barrier(X0s, ys, t0s, t0, t0 + 1.0, level, "0");
}

/* edges X0 y0 t0 N1 Stail : upper edge y = 1 (N0 <= N <= N1), tail region N >= N1 and right edge N = N1.
   Stail is the arb upper bound A + max_N B - 1 for |f_t - 1| on N >= N1 (full_run.part_tail). */
static int run_edges(const char *X0s, const char *ys, const char *ts, ulong N1, const char *Stails)
{
    arb_init(PI); arb_const_pi(PI, P);
    arb_t X0, x, y, t, one, tmp, tmp2, lN0, Ehi, Elo, Ear, need, eb, S, St, sig1, ln, w, s1, s2, s3, ph, rr, xN1, eN1;
    arb_init(X0); arb_init(x); arb_init(y); arb_init(t); arb_init(one); arb_init(tmp); arb_init(tmp2); arb_init(lN0);
    arb_init(Ehi); arb_init(Elo); arb_init(Ear); arb_init(need); arb_init(eb); arb_init(S); arb_init(St); arb_init(sig1);
    arb_init(ln); arb_init(w); arb_init(s1); arb_init(s2); arb_init(s3); arb_init(ph); arb_init(rr); arb_init(xN1); arb_init(eN1);
    arb_set_str(X0, X0s, P); arb_set_str(y, ys, P); arb_set_str(t, ts, P); arb_one(one);
    arb_set_d(tmp, 0.5); arb_sub(x, X0, tmp, P);
    arb_div(tmp, x, PI, P); arb_mul_2exp_si(tmp, tmp, -2); arb_div_ui(tmp2, t, 16, P); arb_add(tmp, tmp, tmp2, P);
    arb_sqrt(tmp, tmp, P); arb_floor(tmp, tmp, P);
    fmpz_t f; fmpz_init(f);
    if (!arb_get_unique_fmpz(f, tmp)) { printf("{\"mode\":\"edges\",\"error\":\"N0\"}\n"); return 1; }
    ulong N0 = fmpz_get_ui(f);
    arb_log_ui(lN0, N0, P);
    region_ii_budget(Ehi, Elo, Ear, need, x, y, t, N0);
    e_bound(eb, x, y, t, 1);
    /* upper edge y = 1: Re s* >= 1 + (t/2) log N0 - 1e-6 =: sig1;  |gamma| n^y <= 1.03 n/N (y = 1).
       |f - 1| <= sum_{2<=n<=N1} b_n n^{-sig1} + 1.03/N0 sum_{n<=N0} b_n n^{1-sig1} + 1.03 sum_{N0<n<=N1} b_n n^{-sig1} =: S */
    arb_mul(sig1, t, lN0, P); arb_mul_2exp_si(sig1, sig1, -1); arb_add_ui(sig1, sig1, 1, P);
    arb_set_str(tmp, "1e-6", P); arb_sub(sig1, sig1, tmp, P);
    arb_zero(s1); arb_zero(s2); arb_zero(s3);
    for (ulong n = 1; n <= N1; n++) {
        arb_log_ui(ln, n, P);
        arb_sqr(tmp, ln, P); arb_mul(tmp, tmp, t, P); arb_mul_2exp_si(tmp, tmp, -2); arb_submul(tmp, sig1, ln, P);
        arb_exp(w, tmp, P);
        if (n >= 2) arb_add(s1, s1, w, P);
        if (n <= N0) arb_addmul_ui(s2, w, n, P); else arb_add(s3, s3, w, P);
    }
    arb_set_str(tmp2, "1.03", P);
    arb_div_ui(tmp, s2, N0, P); arb_mul(tmp, tmp, tmp2, P); arb_add(S, s1, tmp, P);
    arb_mul(tmp, s3, tmp2, P); arb_add(S, S, tmp, P);
    int up_ok = 0; double up_dist = -1;
    if (arb_lt(S, one)) {
        arb_asin(ph, S, P); arb_add(ph, ph, Ear, P);
        arb_sub(rr, one, S, P); arb_mul(rr, rr, Elo, P);
        up_dist = dist_lb(rr, ph); up_ok = up_dist > ub(need);
    }
    /* tail region N >= N1 and right edge: |f| >= 1 - Stail > e(x_{N1}); dist(E f) > E_hi e(x_{N1}) */
    arb_set_str(St, Stails, P);
    arb_set_ui(xN1, N1); arb_sqr(xN1, xN1, P); arb_mul(xN1, xN1, PI, P); arb_mul_2exp_si(xN1, xN1, 2);
    arb_mul(tmp, PI, t, P); arb_mul_2exp_si(tmp, tmp, -2); arb_sub(xN1, xN1, tmp, P);
    e_bound(eN1, xN1, y, t, 0);
    arb_sub(rr, one, St, P);
    int tail_ok = arb_gt(rr, eN1), right_ok = 0; double right_dist = -1;
    if (arb_lt(St, one)) {
        arb_asin(ph, St, P); arb_add(ph, ph, Ear, P); arb_mul(rr, rr, Elo, P);
        right_dist = dist_lb(rr, ph); arb_mul(tmp, Ehi, eN1, P); right_ok = right_dist > ub(tmp);
    }
    printf("{\"mode\":\"edges\",\"N0\":%lu,\"N1\":%lu,\"E_hi\":%.10e,\"E_lo\":%.10e,\"E_arg\":%.10e,\"need\":%.6e,\"e_barrier\":%.10e,"
           "\"upper_S\":%.10e,\"upper_dist\":%.6e,\"upper_ok\":%s,\"tail_S\":%.10e,\"e_N1\":%.3e,\"tail_ok\":%s,\"right_dist\":%.6e,\"right_ok\":%s}\n",
           N0, N1, ub(Ehi), lb(Elo), ub(Ear), ub(need), ub(eb), ub(S), up_dist, up_ok ? "true" : "false",
           ub(St), ub(eN1), tail_ok ? "true" : "false", right_dist, right_ok ? "true" : "false");
    fflush(stdout);
    return (up_ok && tail_ok && right_ok) ? 0 : 2;
}

int main(int argc, char **argv)
{
    setup_primes();
    if (argc >= 7 && !strcmp(argv[1], "edges"))
        return run_edges(argv[2], argv[3], argv[4], strtoul(argv[5], 0, 10), argv[6]);
    if (argc >= 6 && !strcmp(argv[1], "leftedge"))
        return run_leftedge(argv[2], argv[3], argv[4], atoi(argv[5]));
    if (argc >= 6 && !strcmp(argv[1], "block"))
        return run_block(strtoul(argv[2], 0, 10), strtoul(argv[3], 0, 10), argv[4], argv[5],
                         argc >= 8 ? strtoul(argv[6], 0, 10) : 2, argc >= 8 ? strtoul(argv[7], 0, 10) : 0);
    if (argc >= 9 && !strcmp(argv[1], "barrier"))
        return run_barrier(argv[2], argv[3], argv[4], atof(argv[5]), atof(argv[6]), atoi(argv[7]), argv[8]);
    fprintf(stderr, "usage: dbn block Nm Np t0 y0 | dbn barrier X0 y0 t0 ta tb level margin\n");
    return 1;
}

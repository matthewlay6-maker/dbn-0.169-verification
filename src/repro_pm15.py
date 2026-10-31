"""Session 2b: independent reproduction of Polymath15's Lambda <= 0.22 (Section 8 of arXiv:1904.12438).

Parameters: t0 = y0 = 0.2, X0 = 6e10 + 83952, X = X0 - 0.5, N0 = 69098, N1 = 1.5e6.
Claim (i) (RH to height X/2 ~ 3e10) is an external input (Platt 2017; superseded by Platt-Trudgian 2021).
This file reproduces:
  sweep : e_A + e_B + e_{C,0} <= 1.25e-3 on all regions        (arb)
  (c)   : |f_t - 1| <= 0.955 for N >= N1, t = 0.2               (arb for A; B monotone check)
  (b)   : F_{N-,N+} >= 2.15e-3 on intervals covering [N0, N1]   (float64 + explicit rounding allowance)
The barrier claim (a) lives in barrier_pm15.py.
Usage: python repro_pm15.py [sweep|c|b|all]
"""
import json, math, sys, time
import numpy as np
import flint
from flint import arb, acb

flint.ctx.prec = 128
T0 = arb("0.2")
X0 = arb("60000083952")
XL = X0 - arb("0.5")
N0, N1 = 69098, 1_500_000
PI = arb.pi()
OUT = "results/repro_pm15.json"


def load():
    try:
        return json.load(open(OUT))
    except FileNotFoundError:
        return {}


def save(d):
    json.dump(d, open(OUT, "w"), indent=1)


def sweep():
    """Proposition 'sweep' (Section 8): reproduce delta_1, e_A+e_B and e_{C,0} bounds."""
    lx = (XL / (4 * PI)).log()
    d1 = (T0 ** 2 / 16 * lx ** 2 + arb("0.626")) / (XL - arb("6.66"))
    # e_A + e_B <= (e^{d1}-1)(1 + 1.021) F,  F = sum_{n<=N} n^{-0.5999} <= N^{0.4001}/0.4001,
    # N <= sqrt(x/4pi + t/16); d1 * F ~ log^2 x / x^{0.8} is decreasing, so evaluate at x = X.
    # (The paper writes 1.023, 0.5999 and 0.80005 in this chain; recomputed here from its own steps.)
    Nmax = (XL / (4 * PI) + T0 / 16).sqrt()
    eab = (d1.exp() - 1) * (1 + (arb("0.02")).exp() * arb("1.0000001")) * Nmax ** arb("0.4001") / arb("0.4001")
    ymax = arb("0.6").sqrt()
    ec0 = (-(1 + T0) / 4 * lx).exp() * ((3 * abs(acb(lx, PI / 2)) + arb("3.58")) / (XL - arb("8.52"))).exp() \
        * (1 + arb("1.24") * (arb(3) ** ymax + arb(3) ** (-ymax)) / (N0 - arb("0.125")) + arb("6.92") / (XL - arb("6.66")))
    # N range check: N is constant (= N0) on the barrier
    nlo = int(((XL / (4 * PI)).sqrt()).floor().unique_fmpz())
    nhi = int((((X0 + arb("0.5")) / (4 * PI) + T0 / 16).sqrt()).floor().unique_fmpz())
    res = dict(delta1=float(d1.upper()), paper_delta1=3.12e-11,
               eAB=float(eab.upper()), paper_eAB=3.444e-9,
                              eC0=float(ec0.upper()), paper_eC0=1.249e-3,
               total=float((eab + ec0).upper()), claim_total=1.25e-3,
               ok=bool((eab + ec0) < arb("1.25e-3")),
               N_on_barrier=[nlo, nhi])
    print(f"sweep: eAB<={res['eAB']:.4e} eC0<={res['eC0']:.6e} total<={res['total']:.6e} "
          f"(paper 1.25e-3) ok={res['ok']} N on barrier {nlo}..{nhi}")
    return res


def bnt(n, t):
    ln = arb(n).log()
    return (t / 4 * ln * ln).exp()


def claim_c():
    """A at N = N1 (sums to N0) and B(N) on [N1, 10^40] (log grid), t = 0.2."""
    t = T0
    lN = arb(N1).log()
    xN1 = 4 * PI * arb(N1) ** 2 - PI * t / 4
    s1 = arb("0.6") + arb("0.1") * lN - arb("5.1e-29")
    s2 = arb("0.4") + arb("0.1") * lN - arb("7.09e-16")
    A1 = arb(0)
    A2 = arb(0)
    for n in range(1, N0 + 1):
        ln = arb(n).log()
        b = (t / 4 * ln * ln).exp()
        A1 += b * (-s1 * ln).exp()
        A2 += b * (-s2 * ln).exp()
    A = A1 + arb("1.006") * (xN1 / (4 * PI)) ** arb("-0.1") * A2

    def Bof(N):
        N = arb(N)
        lN = N.log()
        s1 = arb("0.6") + arb("0.1") * lN - arb("5.1e-29")
        s2 = arb("0.4") + arb("0.1") * lN - arb("7.09e-16")
        bN0 = (t / 4 * arb(N0).log() ** 2).exp()
        bN = (t / 4 * lN ** 2).exp()
        m1 = max(arb(N0) ** (1 - s1) * bN0, N ** (1 - s1) * bN, key=lambda v: float(v))
        m2 = max(arb(N0) ** (1 - s2) * bN0, N ** (1 - s2) * bN, key=lambda v: float(v))
        return (m1 + arb("1.006") * N ** arb("-0.2") * m2) * (N / N0).log()

    grid = [N1 * 10 ** (k / 10) for k in range(0, 341)]
    Bs = [float(Bof(arb(int(g))).upper()) for g in grid]
    mono = all(Bs[i + 1] <= Bs[i] for i in range(len(Bs) - 1))
    res = dict(A=float(A.upper()), paper_A=1.88, B_at_N1=Bs[0], paper_B=0.075,
               B_max_on_grid=max(Bs), B_decreasing_on_grid=mono,
               total=float(A.upper()) + max(Bs), need_lt=1.955,
               ok=bool(float(A.upper()) + max(Bs) < 1.955))
    print(f"claim (c): A<={res['A']:.4f} (paper 1.88), B(N1)<={res['B_at_N1']:.4f} (paper 0.075), "
          f"B decreasing on grid to 1e40: {mono}, A+B={res['total']:.4f} < 1.955: {res['ok']}")
    return res


def coeffs(Nm, Np, t=0.2):
    """Mollified coefficients for claim (b), rigorous in N in [Nm, Np]:
    beta_n = sum_{d|n, d|30, n/d <= N} lambda_d b_{n/d};  alpha_n = 1.03 sum lambda_d (n/d)^0.2 b_{n/d}.
    Contributions with n/d <= Nm are fixed ('core'); those with Nm < n/d <= Np depend on N and
    are kept as absolute values ('var') and added in worst case."""
    D = 30
    divs = [1, 2, 3, 5, 6, 10, 15, 30]
    lam = {}
    for d in divs:
        v = 1.0
        for p in (2, 3, 5):
            if d % p == 0:
                v *= -math.exp(t / 4 * math.log(p) ** 2)
        lam[d] = v
    L = D * Np
    bc = np.zeros(L + 1)
    ac = np.zeros(L + 1)
    bv = np.zeros(L + 1)
    av = np.zeros(L + 1)
    m = np.arange(1, Np + 1, dtype=np.float64)
    lm = np.log(m)
    b = np.exp(t / 4 * lm * lm)
    a = 1.03 * np.exp(0.2 * lm) * b
    core = m <= Nm
    for d in divs:
        idx = d * np.arange(1, Np + 1)
        bc[idx[core]] += lam[d] * b[core]
        ac[idx[core]] += lam[d] * a[core]
        bv[idx[~core]] += abs(lam[d]) * b[~core]
        av[idx[~core]] += abs(lam[d]) * a[~core]
    return bc, ac, bv, av


def F_interval(Nm, Np):
    bc, ac, bv, av = coeffs(Nm, Np)
    sig = 0.599 + 0.1 * math.log(Nm)
    c = Nm ** -0.2
    a1 = ac[1]
    ratio = (1 - c * a1) / (1 + c * a1)
    n = np.arange(len(bc), dtype=np.float64)
    tot = 0.0
    CH = 5_000_000
    for s in range(2, len(bc), CH):
        e = min(len(bc), s + CH)
        w = np.exp(-sig * np.log(n[s:e]))
        num = np.maximum(np.abs(bc[s:e] - c * ac[s:e]), ratio * np.abs(bc[s:e] + c * ac[s:e]))
        num += bv[s:e] + c * av[s:e]  # worst case for the N-dependent contributions
        tot += float(np.sum(num * w))
    F = 1 - c * a1 - tot
    # float64 allowance: each term has relative error < 1e-13 (exp/log, ~10 flops), the sum
    # of len(bc) positive terms adds < len * 2^-53 relative; bound total error by 1e-6 absolute.
    allow = 1e-6
    return F, F - allow, len(bc)


def claim_b():
    ivs = [(69098, 80000), (80000, 110000), (110000, 220000), (220000, 1500000)]
    paper = [0.0263, 0.0470, 0.093, 0.060]
    rows = []
    for (a, b), p in zip(ivs, paper):
        t0 = time.time()
        F, Flo, L = F_interval(a, b)
        rows.append(dict(Nminus=a, Nplus=b, F=F, F_lower=Flo, paper=p, need=2.15e-3,
                         ok=bool(Flo >= 2.15e-3), terms=L, s=time.time() - t0))
        print(f"claim (b): F[{a},{b}] = {F:.4f} (paper {p}) ok={Flo >= 2.15e-3} {time.time()-t0:.0f}s", flush=True)
    return dict(rows=rows, ok=all(r["ok"] for r in rows))


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    d = load()
    if which in ("sweep", "all"):
        d["sweep"] = sweep(); save(d)
    if which in ("c", "all"):
        d["claim_c"] = claim_c(); save(d)
    if which in ("b", "all"):
        d["claim_b"] = claim_b(); save(d)


def F_block_exact(Nm, Np, t=0.2):
    """min over N in [Nm, Np] of the paper's F, but with the TRUE truncated coefficients
    beta_n^(N) = sum_{d|n, d|30, n/d <= N} lambda_d b_{n/d} (and alpha likewise), updated
    incrementally as N grows (only n = d(N+1), d | 30, change).  sigma_{Nm}, c = Nm^-0.2 fixed."""
    divs = [1, 2, 3, 5, 6, 10, 15, 30]
    lam = {}
    for d in divs:
        v = 1.0
        for p in (2, 3, 5):
            if d % p == 0:
                v *= -math.exp(t / 4 * math.log(p) ** 2)
        lam[d] = v
    bc, ac, _, _ = coeffs(Nm, Nm)
    L = 30 * Np + 1
    beta = np.zeros(L); alph = np.zeros(L)
    beta[:len(bc)] = bc; alph[:len(ac)] = ac
    sig = 0.599 + 0.1 * math.log(Nm)
    c = Nm ** -0.2
    a1 = alph[1]
    r = (1 - c * a1) / (1 + c * a1)
    n = np.arange(L, dtype=np.float64)
    w = np.zeros(L); w[2:] = np.exp(-sig * np.log(n[2:]))
    g = np.maximum(np.abs(beta - c * alph), r * np.abs(beta + c * alph)) * w
    g[:2] = 0
    T = float(np.sum(g))
    Fmin, Nat = 1 - c * a1 - T, Nm
    for m in range(Nm + 1, Np + 1):
        lm = math.log(m)
        bm = math.exp(t / 4 * lm * lm)
        am = 1.03 * math.exp(0.2 * lm) * bm
        for d in divs:
            k = d * m
            T -= g[k]
            beta[k] += lam[d] * bm
            alph[k] += lam[d] * am
            gk = max(abs(beta[k] - c * alph[k]), r * abs(beta[k] + c * alph[k])) * w[k]
            g[k] = gk
            T += gk
        F = 1 - c * a1 - T
        if F < Fmin:
            Fmin, Nat = F, m
    drift = abs(T - float(np.sum(g)))
    return Fmin, Nat, drift, sig


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "bexact":
    blocks = [int(v) for v in sys.argv[2].split(",")] if len(sys.argv) > 2 else [69098, 80000, 110000, 220000, 1500000]
    d = load()
    rows = []
    for a, b in zip(blocks[:-1], blocks[1:]):
        t0 = time.time()
        Fmin, Nat, drift, sig = F_block_exact(a, b)
        ok = Fmin - 1e-6 - drift >= 2.15e-3
        rows.append(dict(Nminus=a, Nplus=b, F_min=Fmin, argmin_N=Nat, drift=drift, ok=bool(ok), s=time.time() - t0))
        print(f"(b) exact: [{a},{b}] min_N F = {Fmin:.5f} at N={Nat} ok={ok} drift={drift:.1e} {time.time()-t0:.0f}s", flush=True)
    d.setdefault("claim_b_exact", {})[",".join(map(str, blocks))] = dict(rows=rows, ok=all(r["ok"] for r in rows))
    save(d)

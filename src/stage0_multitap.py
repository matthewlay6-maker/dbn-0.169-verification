"""Stage 0 (float64 scouting, not rigorous): multi-tap mollifiers.

Mollifier E(s) = prod_p E_p(s),  E_p = 1 + sum_{k=1..K_p} m_{p,k} p^{-k s}   (Euler single tap: K_p = 1, m_{p,1} = -b_p).
Mollified coefficients:  beta_n = sum_{d | n, d in Div} lambda_d b_{n/d} [n/d <= N],  lambda_d = prod_p m_{p, e_p(d)} (m_{p,0} = 1),
and alpha_n likewise with gcon (n/d)^y0.  Bound (Polymath15 section 8 improved triangle inequality, tight constants):
    F(m) = 1 - c*alpha_1 - sum_{n>=2} max(|beta_n - c alpha_n|, r |beta_n + c alpha_n|) n^{-sigma}
Any real m gives a rigorous mollifier (argument principle on E*H/B); the price is max|E| in the error budget.
Sparse evaluation: all (d, m) pairs are generated once, sorted by n = d*m; each F(m) evaluation re-weights them.
"""
import itertools, json, math, sys, time
import numpy as np
from scipy.optimize import minimize

X = 6000347387366 - 0.5


class Setup:
    def __init__(self, N, t0, y0, taps):            # taps: {p: K_p}
        self.N, self.t0, self.y0, self.taps = N, t0, y0, dict(taps)
        self.P = sorted(taps)
        self.sig = (1 + y0) / 2 + t0 / 2 * math.log(N) - 1e-9
        self.g = math.exp(0.02 * y0) * (1 + 1e-9)
        self.c = N ** -y0
        exps = list(itertools.product(*[range(taps[p] + 1) for p in self.P]))
        self.exps = np.array(exps)                       # divisor exponent vectors
        self.dvals = np.array([math.prod(p ** e for p, e in zip(self.P, ev)) for ev in exps], dtype=np.int64)
        m = np.arange(1, N + 1, dtype=np.int64)
        lm = np.log(m.astype(np.float64))
        b = np.exp(t0 / 4 * lm * lm)
        a = self.g * np.exp(y0 * lm) * b
        nn, di, bv, av = [], [], [], []
        for k, d in enumerate(self.dvals):
            nn.append(m * d); di.append(np.full(N, k, np.int32)); bv.append(b); av.append(a)
        nn = np.concatenate(nn); di = np.concatenate(di); bv = np.concatenate(bv); av = np.concatenate(av)
        order = np.argsort(nn, kind="stable")
        self.n = nn[order]; self.di = di[order]; self.bv = bv[order]; self.av = av[order]
        self.starts = np.flatnonzero(np.r_[True, self.n[1:] != self.n[:-1]])
        un = self.n[self.starts]
        self.w = np.exp(-self.sig * np.log(un.astype(np.float64)))
        self.is1 = (un == 1)

    def default_m(self):
        """single Euler taps: m_{p,1} = -b_p, higher taps 0."""
        out = {}
        for p in self.P:
            out[p] = [-math.exp(self.t0 / 4 * math.log(p) ** 2)] + [0.0] * (self.taps[p] - 1)
        return out

    def lam(self, mt):
        lam = np.ones(len(self.dvals))
        for j, p in enumerate(self.P):
            coeff = np.array([1.0] + list(mt[p]))
            lam *= coeff[self.exps[:, j]]
        return lam

    def F(self, mt):
        lam = self.lam(mt)
        beta = np.add.reduceat(lam[self.di] * self.bv, self.starts)
        alph = np.add.reduceat(lam[self.di] * self.av, self.starts)
        c = self.c
        a1 = alph[self.is1][0]
        r = (1 - c * a1) / (1 + c * a1)
        term = np.maximum(np.abs(beta - c * alph), r * np.abs(beta + c * alph)) * self.w
        term[self.is1] = 0.0
        return 1 - c * a1 - float(np.sum(term))

    def Ehi(self, mt):
        v = 1.0
        for p in self.P:
            z = p ** -self.sig
            v *= 1 + sum(abs(mk) * z ** (k + 1) for k, mk in enumerate(mt[p]))
        return v


def optimise(S, x0_m, maxiter=600):
    keys = [(p, k) for p in S.P for k in range(S.taps[p])]
    x0 = np.array([x0_m[p][k] for p, k in keys])

    def unpack(x):
        mt = {p: [0.0] * S.taps[p] for p in S.P}
        for (p, k), v in zip(keys, x):
            mt[p][k] = v
        return mt
    f = lambda x: -S.F(unpack(x))
    r = minimize(f, x0, method="Nelder-Mead", options=dict(maxiter=maxiter, xatol=1e-6, fatol=1e-8, adaptive=True))
    return -r.fun, unpack(r.x), r.nfev


if __name__ == "__main__":
    lam_t, y0 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.170, 0.08
    t0 = lam_t - y0 * y0 / 2
    N0 = int(math.floor(math.sqrt(X / (4 * math.pi) + t0 / 16)))
    out = dict(Lambda=lam_t, y0=y0, t0=t0, N0=N0, variants={}, halving={})
    T = time.time()
    variants = {"A single exact (Euler-11)": {2: 1, 3: 1, 5: 1, 7: 1, 11: 1},
                "B single optimised": {2: 1, 3: 1, 5: 1, 7: 1, 11: 1},
                "C 2-tap on 2,3": {2: 2, 3: 2, 5: 1, 7: 1, 11: 1},
                "D 3-tap on 2, 2-tap on 3": {2: 3, 3: 2, 5: 1, 7: 1, 11: 1}}
    for name, taps in variants.items():
        S = Setup(N0, t0, y0, taps)
        m0 = S.default_m()
        F0 = S.F(m0)
        if name.startswith("A"):
            F, m, nf = F0, m0, 0
        else:
            F, m, nf = optimise(S, m0)
        out["variants"][name] = dict(F=F, F_start=F0, taps={str(p): m[p] for p in S.P}, Ehi=S.Ehi(m), evals=nf)
        print(f"{name:28s}: F = {F:+.5f} (start {F0:+.5f}, max|E| = {S.Ehi(m):.3f}, {nf} evals, {time.time()-T:.0f}s)", flush=True)
    # halving law: add primes one at a time, single exact taps vs optimised 2-taps (3-tap on 2)
    primes = [2, 3, 5, 7, 11]
    for kind in ("single exact", "multi-tap optimised"):
        seq = []
        for k in range(1, len(primes) + 1):
            taps = {p: (1 if kind == "single exact" else (3 if p == 2 else 2)) for p in primes[:k]}
            S = Setup(N0, t0, y0, taps)
            m0 = S.default_m()
            F = S.F(m0) if kind == "single exact" else optimise(S, m0, maxiter=400)[0]
            seq.append(F)
            print(f"halving [{kind}] primes<= {primes[k-1]:2d}: F = {F:+.5f}" + (f"  gain {F - seq[-2]:+.5f}" if k > 1 else ""), flush=True)
        out["halving"][kind] = seq
    json.dump(out, open(f"results/stage0_multitap_{lam_t}.json", "w"), indent=1)

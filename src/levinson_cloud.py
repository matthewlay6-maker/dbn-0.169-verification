"""Stage 0b (float64 scouting, not rigorous): Levinson/Selberg-type long mollifier for the lower-edge bound.

M(s) = sum_{d <= y, squarefree} lambda_d d^{-s},  lambda_d = mu(d) P(log(y/d)/log y),  P(0) = 0, P(1) = 1.
Mollified coefficients at N = N0 (single N, exact truncation):
  beta_n = sum_{d | n, d <= y} lambda_d b_{n/d} [n/d <= N0],   alpha_n = sum lambda_d gcon (n/d)^{y0} b_{n/d} [...]
  F = 1 - c alpha_1 - sum_{n>=2} max(|beta_n - c alpha_n|, r |beta_n + c alpha_n|) n^{-sigma}     (tight constants)
Also reports max|M| <= sum |lambda_d| d^{-sigma} and the lower bound 1 - sum_{d>=2} |lambda_d| d^{-sigma} for |M|.
Task = CLOUD_RUN_TASK_INDEX into /app/specs.json; writes /out/levinson/task-XXX.json.
"""
import itertools, json, math, os, time
import numpy as np
from scipy.optimize import minimize

X = 6000347387366 - 0.5


def mobius_upto(y):
    mu = np.ones(y + 1, dtype=np.int64); mu[0] = 0
    isp = np.ones(y + 1, bool); isp[:2] = False
    for p in range(2, y + 1):
        if isp[p]:
            isp[2 * p::p] = False
            mu[p::p] *= -1
            mu[p * p::p * p] = 0
    return mu


class Lev:
    def __init__(self, N, t0, y0, ymax):
        self.N, self.t0, self.y0, self.ymax = N, t0, y0, ymax
        self.sig = (1 + y0) / 2 + t0 / 2 * math.log(N) - 1e-9
        self.g = math.exp(0.02 * y0) * (1 + 1e-9)
        self.c = N ** -y0
        mu = mobius_upto(ymax)
        self.ds = np.array([d for d in range(1, ymax + 1) if mu[d] != 0])
        self.mu = mu[self.ds]
        m = np.arange(1, N + 1, dtype=np.float64)
        lm = np.log(m)
        self.b = np.exp(t0 / 4 * lm * lm)
        self.a = self.g * np.exp(y0 * lm) * self.b

    def bw(self):
        """b-weights prod_{p | d} b_p (Euler-like), so the mollifier cancels the b-weighted first echoes."""
        out = np.ones(len(self.ds))
        for i, d in enumerate(self.ds):
            v, q = int(d), 2
            while v > 1:
                if v % q == 0:
                    out[i] *= math.exp(self.t0 / 4 * math.log(q) ** 2); v //= q
                q += 1
        return out

    def lam(self, coeffs, taper="poly", weighted=True):
        """lambda_d = mu(d) [prod_{p|d} b_p] P(x), x = log(y/d)/log y.
        taper 'poly': P(x) = x + a2 (x^2 - x) + a3 (x^3 - x);  'flat': P = 1 for d < y (sharp cut);
        'step': P(x) = min(1, x / a2) (flat then linear to 0 near d = y)."""
        a2, a3 = coeffs
        x = np.log(self.ymax / self.ds) / math.log(self.ymax)
        if taper == "flat":
            P = np.ones(len(self.ds))
        elif taper == "step":
            P = np.minimum(1.0, x / max(a2, 1e-6))
        else:
            P = x + a2 * (x * x - x) + a3 * (x ** 3 - x)
        lam = self.mu * P * (self.bw() if weighted else 1.0)
        lam[self.ds == 1] = 1.0
        return lam

    def F(self, lam):
        N, c, g = self.N, self.c, self.g
        r = (1 - c * g) / (1 + c * g)
        L = self.ymax * N
        tot = 0.0
        CH = 16_000_000
        for s in range(2, L + 1, CH):
            e = min(L + 1, s + CH)
            bc = np.zeros(e - s); ac = np.zeros(e - s)
            for d, ld in zip(self.ds, lam):
                if ld == 0.0:
                    continue
                mlo = -(-s // d); mhi = min((e - 1) // d, N)
                if mhi < mlo:
                    continue
                idx = np.arange(mlo, mhi + 1) * d - s
                bc[idx] += ld * self.b[mlo - 1:mhi]
                ac[idx] += ld * self.a[mlo - 1:mhi]
            n = np.arange(s, e, dtype=np.float64)
            w = np.exp(-self.sig * np.log(n))
            tot += float(np.sum(np.maximum(np.abs(bc - c * ac), r * np.abs(bc + c * ac)) * w))
        return 1 - c * g - tot

    def Mbounds(self, lam):
        w = self.ds.astype(np.float64) ** -self.sig
        hi = float(np.sum(np.abs(lam) * w))
        lo = 1.0 - float(np.sum((np.abs(lam) * w)[self.ds > 1]))
        return hi, lo


def euler11_F(N, t0, y0):
    P = (2, 3, 5, 7, 11)
    D = math.prod(P)
    divs = sorted(math.prod(cc) for k in range(len(P) + 1) for cc in itertools.combinations(P, k))
    lam = {d: math.prod(-math.exp(t0 / 4 * math.log(p) ** 2) for p in P if d % p == 0) for d in divs}
    S = Lev.__new__(Lev); S.N, S.t0, S.y0, S.ymax = N, t0, y0, D
    S.sig = (1 + y0) / 2 + t0 / 2 * math.log(N) - 1e-9
    S.g = math.exp(0.02 * y0) * (1 + 1e-9); S.c = N ** -y0
    m = np.arange(1, N + 1, dtype=np.float64); lm = np.log(m)
    S.b = np.exp(t0 / 4 * lm * lm); S.a = S.g * np.exp(y0 * lm) * S.b
    S.ds = np.array(divs); S.mu = np.ones(len(divs))
    return S.F(np.array([lam[d] for d in divs]))


if __name__ == "__main__":
    idx = int(os.environ.get("CLOUD_RUN_TASK_INDEX", "0"))
    sp = json.load(open("/app/specs.json"))[idx]
    lam_t, y0, ymax, mode = sp["Lambda"], sp["y0"], sp["y"], sp["mode"]
    t0 = lam_t - y0 * y0 / 2
    N0 = int(math.floor(math.sqrt(X / (4 * math.pi) + t0 / 16)))
    T = time.time()
    if mode == "euler11":
        res = dict(F=euler11_F(N0, t0, y0), evals=1)
    else:
        S = Lev(N0, t0, y0, ymax)
        if mode == "flat":
            best = (0.0, 0.0); lam = S.lam(best, "flat"); F = S.F(lam); nf = 1; taper = "flat"
        elif mode == "step_opt":
            from scipy.optimize import minimize_scalar
            r = minimize_scalar(lambda a: -S.F(S.lam((a, 0.0), "step")), bounds=(0.02, 1.0), method="bounded",
                                options=dict(maxiter=sp.get("maxiter", 25), xatol=1e-3))
            best = (float(r.x), 0.0); F = -r.fun; nf = r.nfev; lam = S.lam(best, "step"); taper = "step"
        else:
            f = lambda x: -S.F(S.lam(x, "poly"))
            r = minimize(f, np.zeros(2), method="Nelder-Mead", options=dict(maxiter=sp.get("maxiter", 60), xatol=1e-4, fatol=1e-7))
            best = tuple(r.x); F = -r.fun; nf = r.nfev; lam = S.lam(best, "poly"); taper = "poly"
        hi, lo = S.Mbounds(lam)
        res = dict(F=F, taper=taper, P_coeffs=list(best), M_hi=hi, M_lo_bound=lo, n_terms=int(len(S.ds)), evals=nf)
    res.update(spec=sp, N0=N0, t0=t0, seconds=time.time() - T)
    os.makedirs("/out/levinson", exist_ok=True)
    json.dump(res, open(f"/out/levinson/task-{idx:03d}.json", "w"), indent=1)
    print(json.dumps(res), flush=True)

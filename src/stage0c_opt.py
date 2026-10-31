"""Stage 0c Program A (float64 scouting, not rigorous): globally optimal mollifier weights for the triangle floor.

For a fixed support S (with 1 in S, mu_1 = 1) the hardest-window floor at N = N0 (single N, improved triangle weight)
    F(mu) = 1 - c g - sum_{n>=2} w_n max(|u_n|, |v_n|),   u = beta - c alpha,  v = r (beta + c alpha),
    beta_n = sum_{m in S, m|n, n/m <= N} mu_m b_{n/m},   alpha_n = sum mu_m a_{n/m},
is concave in mu (beta, alpha linear), so maximising it is convex.  Solved by L-BFGS on a smoothed objective
(|z| -> sqrt(z^2 + e^2), max -> (p + q + sqrt((p-q)^2 + e^2))/2 on the weighted scale) with e -> 0, exact gradients.
The reported F is always the exact unsmoothed value.  Also reports the error charge max|M| * e_total.
Task = CLOUD_RUN_TASK_INDEX into /app/specs.json; writes /out/stage0c/task-XXX.json.
"""
import json, math, os, time
import numpy as np
from scipy.optimize import minimize

X = 6000347387366 - 0.5


def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s)


def factor(m):
    out, q = {}, 2
    while m > 1:
        while m % q == 0:
            out[q] = out.get(q, 0) + 1; m //= q
        q += 1
    return out


def support(kind, size):
    """(i) all m <= B; (ii) 13-smooth; (iii) 17-smooth with exponents > 1 only for 2, 3; (iv) squarefree <= 300."""
    if kind == "control":
        return np.array([m for m in range(1, 301) if all(e == 1 for e in factor(m).values())], dtype=np.int64)
    out, m = [], 1
    while len(out) < size:
        f = factor(m)
        ok = True
        if kind == "smooth13":
            ok = all(p <= 13 for p in f)
        elif kind == "smooth17":
            ok = all(p <= 17 for p in f) and all(e == 1 for p, e in f.items() if p > 3)
        if ok:
            out.append(m)
        m += 1
    return np.array(out, dtype=np.int64)


def init_mu(ms, t0, kind, coeffs=(-1.805, 0.0), y=300):
    """b-weighted Moebius start (exactly the heuristic long mollifier for the control)."""
    mu = np.zeros(len(ms))
    for i, m in enumerate(ms):
        f = factor(int(m))
        if any(e > 1 for e in f.values()):
            continue
        bw = math.prod(math.exp(t0 / 4 * math.log(p) ** 2) for p in f)
        if kind == "control":
            x = math.log(y / m) / math.log(y)
            P = x + coeffs[0] * (x * x - x) + coeffs[1] * (x ** 3 - x)
        else:
            x = math.log(ms.max() / m) / math.log(ms.max()) if ms.max() > 1 else 1.0
            P = x + coeffs[0] * (x * x - x)
        mu[i] = (-1) ** len(f) * bw * P
    mu[ms == 1] = 1.0
    return mu


class Floor:
    def __init__(self, N, t0, y0, ms):
        self.N, self.t0, self.y0, self.ms = N, t0, y0, ms
        self.sig = (1 + y0) / 2 + t0 / 2 * math.log(N) - 1e-9
        self.g = math.exp(0.02 * y0) * (1 + 1e-9)
        self.c = N ** -y0
        self.r = (1 - self.c * self.g) / (1 + self.c * self.g)
        k = np.arange(1, N + 1, dtype=np.float64); lk = np.log(k)
        self.b = np.exp(t0 / 4 * lk * lk)
        self.a = self.g * np.exp(y0 * lk) * self.b
        self.bm = self.b - self.c * self.a          # d u / d mu_m pattern
        self.bp = self.r * (self.b + self.c * self.a)  # d v / d mu_m pattern
        self.L = int(ms.max()) * N

    def chunk_starts(self, CH=16_000_000):
        return list(range(2, self.L + 1, CH))

    def value_grad(self, mu, eps, starts=None, CH=16_000_000, raw=False):
        """returns (F_exact, F_smooth, grad_smooth wrt mu); with raw=True the partial sums over the given chunks."""
        N, ms = self.N, self.ms
        grad = np.zeros(len(ms)); tot_exact = 0.0; tot_s = 0.0
        for s in (self.chunk_starts(CH) if starts is None else starts):
            e = min(self.L + 1, s + CH)
            u = np.zeros(e - s); v = np.zeros(e - s)
            sl = []
            for j, (m, mj) in enumerate(zip(ms, mu)):
                klo = -(-s // m); khi = min((e - 1) // m, N)
                if khi < klo:
                    sl.append(None); continue
                idx = np.arange(klo, khi + 1) * m - s
                sl.append((idx, klo, khi))
                if mj != 0.0:
                    u[idx] += mj * self.bm[klo - 1:khi]
                    v[idx] += mj * self.bp[klo - 1:khi]
            w = np.exp(-self.sig * np.log(np.arange(s, e, dtype=np.float64)))
            U, V = w * u, w * v
            tot_exact += float(np.sum(np.maximum(np.abs(U), np.abs(V))))
            aU = np.sqrt(U * U + eps * eps); aV = np.sqrt(V * V + eps * eps)
            dlt = np.sqrt((aU - aV) ** 2 + eps * eps)
            tot_s += float(np.sum(0.5 * (aU + aV + dlt)))
            # d/dU of smoothed max: 0.5 (1 + (aU-aV)/dlt) * U/aU ; similarly V
            gU = 0.5 * (1 + (aU - aV) / dlt) * U / aU * w
            gV = 0.5 * (1 - (aU - aV) / dlt) * V / aV * w
            for j, it in enumerate(sl):
                if it is None:
                    continue
                idx, klo, khi = it
                grad[j] += float(np.dot(gU[idx], self.bm[klo - 1:khi]) + np.dot(gV[idx], self.bp[klo - 1:khi]))
        if raw:
            return tot_exact, tot_s, grad
        base = 1 - self.c * self.g
        return base - tot_exact, base - tot_s, -grad

    def Mhi(self, mu):
        return float(np.sum(np.abs(mu) * self.ms.astype(np.float64) ** -self.sig))


_FL = None


def _part(args):
    mu, eps, starts = args
    return _FL.value_grad(mu, eps, starts=starts, raw=True)


class ParFloor:
    """Splits the chunk sum of Floor.value_grad over worker processes (fork; workers inherit the Floor)."""
    def __init__(self, fl, nproc):
        global _FL
        import multiprocessing as mp
        _FL = fl
        self.fl = fl
        st = fl.chunk_starts()
        self.nproc = max(1, min(nproc, len(st)))
        self.groups = [st[i::self.nproc] for i in range(self.nproc)]
        self.pool = mp.get_context("fork").Pool(self.nproc) if self.nproc > 1 else None
        self.ms = fl.ms

    def value_grad(self, mu, eps):
        if self.pool is None:
            return self.fl.value_grad(mu, eps)
        parts = self.pool.map(_part, [(mu, eps, g) for g in self.groups])
        te = sum(p[0] for p in parts); ts = sum(p[1] for p in parts); g = sum(p[2] for p in parts)
        base = 1 - self.fl.c * self.fl.g
        return base - te, base - ts, -g

    def Mhi(self, mu):
        return self.fl.Mhi(mu)


def e_total(N, t0, y0):
    x = 4 * math.pi * N * N - math.pi * t0 / 4
    lx = math.log(x / (4 * math.pi))
    d1 = (t0 * t0 / 16 * lx * lx + 0.626) / (x - 6.66)
    s = (1 - y0) / 2
    eab = (math.exp(d1) - 1) * (1 + math.exp(0.02)) * (N + 1) ** s / s
    ec0 = math.exp(-(1 + y0) / 4 * lx) * math.exp(-t0 / 16 * lx * lx + (3 * abs(complex(lx, math.pi / 2)) + 3.58) / (x - 8.52)) \
        * (1 + 1.24 * (3 ** y0 + 3 ** -y0) / (N - 0.125) + 6.92 / (x - 6.66))
    return eab + ec0


if __name__ == "__main__":
    idx = int(os.environ.get("CLOUD_RUN_TASK_INDEX", "0"))
    sp = json.load(open("/app/specs.json"))[idx]
    lam_t, y0 = sp["Lambda"], sp["y0"]
    t0 = lam_t - y0 * y0 / 2
    N0 = int(math.floor(math.sqrt(X / (4 * math.pi) + t0 / 16)))
    T = time.time()
    ms = support(sp["support"], sp["size"])
    base_fl = Floor(N0, t0, y0, ms)
    Fl = ParFloor(base_fl, int(os.environ.get("NPROC", "1")))
    mu0 = init_mu(ms, t0, sp["support"])
    F0, _, _ = Fl.value_grad(mu0, 1e-12)
    F0_serial = base_fl.value_grad(mu0, 1e-12)[0]
    if abs(F0 - F0_serial) > 1e-12:
        raise SystemExit(f"parallel/serial mismatch {F0} vs {F0_serial}")
    free = ms != 1
    x = mu0[free].copy()
    hist = []
    for eps in sp.get("eps", [1e-9, 1e-10, 1e-11]):
        def fg(z):
            mu = mu0.copy(); mu[free] = z
            Fe, Fs, g = Fl.value_grad(mu, eps)
            return -Fs, -g[free]
        r = minimize(fg, x, jac=True, method="L-BFGS-B", options=dict(maxiter=sp.get("maxiter", 60)))
        x = r.x
        mu = mu0.copy(); mu[free] = x
        Fe = Fl.value_grad(mu, 1e-12)[0]
        hist.append(dict(eps=eps, F_exact=Fe, nit=r.nit, nfev=r.nfev))
    mu = mu0.copy(); mu[free] = x
    Ffin = Fl.value_grad(mu, 1e-12)[0]
    et = e_total(N0, t0, y0); Mh = Fl.Mhi(mu)
    top = sorted(zip(np.abs(mu), ms, mu), reverse=True)[:10]
    res = dict(spec=sp, t0=t0, N0=N0, n_support=int(len(ms)), max_m=int(ms.max()), F_start=F0, F_opt=Ffin,
               n_nonzero=int(np.sum(np.abs(mu) > 1e-12)), top10=[(int(m), float(v)) for _, m, v in top],
               Mhi=Mh, error_charge=Mh * et, budget_active=bool(Mh * et > 1e-3 * max(abs(Ffin), 1e-12)),
               history=hist, mu={int(m): float(v) for m, v in zip(ms, mu)}, seconds=time.time() - T)
    os.makedirs("/out/stage0c", exist_ok=True)
    json.dump(res, open(f"/out/stage0c/task-{idx:03d}.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ("mu", "history")}), flush=True)

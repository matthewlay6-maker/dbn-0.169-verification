"""Session 3B: cost map for the Polymath15 method.  For RH verified to height X/2, what is the smallest
Lambda = t0 + y0^2/2 for which the Euler-5-mollified, uniform-in-T lower bound (Section 8, claim (b))
certifies |H_t0/B_t0| > 0 at N = N0 = floor(sqrt(X/4pi))?

Proxy, stated plainly: this checks only the hardest block start N = N0 (F grows with N in all checks so
far) and only hypothesis (ii) there.  The barrier (iii) and the large-N tail are not run for these X.
Criterion A: F(N0; t0, y0) >= |E|max * e(X, t0, y0); criterion B adds a 0.02 safety margin
(Polymath15's own verified F at X=6e10 was 0.0202). Bound e, where e = e_A + e_B + e_{C,0} bound (Theorem 1.3)
and 1.635 bounds |E_{t,5}| (recomputed per t0, y0).
"""
import json, math, sys, time
import numpy as np

DIVS = [1, 2, 3, 5, 6, 10, 15, 30]
OUT = "results/cost_map.json"


def F_at(N, t0, y0):
    lam = {}
    for d in DIVS:
        v = 1.0
        for p in (2, 3, 5):
            if d % p == 0:
                v *= -math.exp(t0 / 4 * math.log(p) ** 2)
        lam[d] = v
    L = 30 * N + 1
    beta = np.zeros(L); alph = np.zeros(L)
    m = np.arange(1, N + 1, dtype=np.float64)
    lm = np.log(m)
    b = np.exp(t0 / 4 * lm * lm)
    a = 1.03 * np.exp(y0 * lm) * b
    for d in DIVS:
        beta[d::d][:N] += lam[d] * b
        alph[d::d][:N] += lam[d] * a
    sig = (1 + y0) / 2 - 0.001 + t0 / 2 * math.log(N)
    c = N ** -y0
    a1 = alph[1]
    r = (1 - c * a1) / (1 + c * a1)
    tot = 0.0
    CH = 20_000_000
    for s in range(2, L, CH):
        e = min(L, s + CH)
        n = np.arange(s, e, dtype=np.float64)
        w = np.exp(-sig * np.log(n))
        tot += float(np.sum(np.maximum(np.abs(beta[s:e] - c * alph[s:e]), r * np.abs(beta[s:e] + c * alph[s:e])) * w))
    return 1 - c * a1 - tot, sig


def err_bound(X, t0, y0, N):
    lx = math.log(X / (4 * math.pi))
    d1 = (t0 ** 2 / 16 * lx ** 2 + 0.626) / (X - 6.66)
    eab = (math.exp(d1) - 1) * (1 + math.exp(0.02)) * (N + 1) ** (0.5 - y0 / 2) / (0.5 - y0 / 2)
    ec0 = (X / (4 * math.pi)) ** (-(1 + y0) / 4) * math.exp(
        -t0 / 16 * lx ** 2 + 1.24 * (3 ** y0 + 3 ** -y0) / (N - 0.125) + (3 * abs(complex(lx, math.pi / 2)) + 10.44) / (X - 12))
    return eab + ec0


def mollifier_max(t0, y0, sig):
    v = 1.0
    for p in (2, 3, 5):
        v *= 1 + math.exp(t0 / 4 * math.log(p) ** 2) / p ** sig
    return v


def feasible(X, t0, y0, margin=0.0):
    N = int(math.floor(math.sqrt(X / (4 * math.pi) + t0 / 16)))
    F, sig = F_at(N, t0, y0)
    need = mollifier_max(t0, y0, sig) * err_bound(X, t0, y0, N)
    return F >= need + margin, F, need, N


def best_lambda(X, y_grid, margin=0.0):
    best = None
    rows = []
    for y0 in y_grid:
        lo, hi = 0.02, 0.35
        ok_hi, *_ = feasible(X, hi, y0, margin)
        if not ok_hi:
            rows.append(dict(y0=y0, t0=None))
            continue
        for _ in range(9):
            mid = (lo + hi) / 2
            if feasible(X, mid, y0, margin)[0]:
                hi = mid
            else:
                lo = mid
        lam = hi + y0 * y0 / 2
        _, F, need, N = feasible(X, hi, y0, margin)
        rows.append(dict(y0=y0, t0=hi, Lambda=lam, F=F, need=need, N0=N))
        if best is None or lam < best["Lambda"]:
            best = rows[-1]
    return best, rows


if __name__ == "__main__":
    Xs = [float(v) for v in sys.argv[1].split(",")] if len(sys.argv) > 1 else \
        [6e10, 2e11, 6e11, 2e12, 6e12, 2e13, 6e13, 2e14, 6e14]
    y_grid = [0.04, 0.06, 0.08, 0.10, 0.12, 0.15]
    try:
        res = json.load(open(OUT))
    except FileNotFoundError:
        res = {}
    for X in Xs:
        key = f"{X:.3g}"
        if key in res:
            continue
        t0 = time.time()
        best, rows = best_lambda(X, y_grid)
        bestS, rowsS = best_lambda(X, y_grid, margin=0.02)
        res[key] = dict(X=X, best=best, rows=rows, best_safety=bestS, rows_safety=rowsS,
                        runtime_s=time.time() - t0)
        json.dump(res, open(OUT, "w"), indent=1)
        print(f"X={X:.3g}: best Lambda={best['Lambda']:.4f} (t0={best['t0']:.4f}, y0={best['y0']}) "
              f"| with 0.02 margin: {bestS['Lambda'] if bestS else None} {time.time()-t0:.0f}s", flush=True)

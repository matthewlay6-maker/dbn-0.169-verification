"""Float scouting: hardest-window (N = N0, single N) lower-edge bound under two weights, same mollifier coefficients.
  plain    (Gomila's certified -DTRIANGLE_WEIGHT): F = 1 - c g - sum_{n>=2} (|beta| + c|alpha|) w
  improved (Polymath15 sec. 8, our pipeline)     : F = 1 - c g - sum_{n>=2} max(|beta - c alpha|, r|beta + c alpha|) w
Mollifiers: Euler-11 (calibration) and the b-weighted long mollifier (y = 300, polynomial taper)."""
import json, math, sys, time
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, "levinson")
import levinson_cloud as L

X = L.X
LONG_COEFFS = (-1.805, 0.0)   # best y = 300 taper found at Lambda = 0.170 (results/lev1)


def both(S, lam):
    N, c, g = S.N, S.c, S.g
    r = (1 - c * g) / (1 + c * g)
    Lmax = int(S.ds.max()) * N
    tp = ti = 0.0
    CH = 16_000_000
    for s in range(2, Lmax + 1, CH):
        e = min(Lmax + 1, s + CH)
        bc = np.zeros(e - s); ac = np.zeros(e - s)
        for d, ld in zip(S.ds, lam):
            mlo = -(-s // d); mhi = min((e - 1) // d, N)
            if mhi < mlo or ld == 0.0:
                continue
            idx = np.arange(mlo, mhi + 1) * d - s
            bc[idx] += ld * S.b[mlo - 1:mhi]
            ac[idx] += ld * S.a[mlo - 1:mhi]
        w = np.exp(-S.sig * np.log(np.arange(s, e, dtype=np.float64)))
        tp += float(np.sum((np.abs(bc) + c * np.abs(ac)) * w))
        ti += float(np.sum(np.maximum(np.abs(bc - c * ac), r * np.abs(bc + c * ac)) * w))
    return 1 - c * g - tp, 1 - c * g - ti


def task(args):
    kind, lam_t, y0 = args
    t0 = lam_t - y0 * y0 / 2
    N0 = int(math.floor(math.sqrt(X / (4 * math.pi) + t0 / 16)))
    S = L.Lev(N0, t0, y0, 300)
    if kind == "long":
        lam = S.lam(LONG_COEFFS, "poly")
    else:   # Euler-11 as an explicit divisor list
        import itertools
        P = (2, 3, 5, 7, 11)
        ds = sorted(math.prod(cc) for k in range(6) for cc in itertools.combinations(P, k))
        S.ds = np.array(ds); S.mu = np.ones(len(ds))
        lam = np.array([math.prod(-math.exp(t0 / 4 * math.log(p) ** 2) for p in P if d % p == 0) for d in ds])
    Fp, Fi = both(S, lam)
    return dict(kind=kind, Lambda=lam_t, y0=y0, t0=t0, N0=N0, F_plain=Fp, F_improved=Fi)


if __name__ == "__main__":
    T = time.time()
    jobs = [("euler11", 0.1787854, math.sqrt(87677 / 2500000))]
    jobs += [("long", lt, y0) for lt in (0.169, 0.170, 0.172, 0.174, 0.176, 0.178) for y0 in (0.06, 0.08, 0.10, 0.14, 0.187)]
    with Pool(8) as p:
        res = p.map(task, jobs)
    json.dump(res, open("results/weights_compare.json", "w"), indent=1)
    cal = res[0]
    print(f"calibration Euler-11 at Gomila row: plain F={cal['F_plain']:+.3e} improved F={cal['F_improved']:+.4f}")
    for kind in ("long",):
        for lt in (0.169, 0.170, 0.172, 0.174, 0.176, 0.178):
            rs = [r for r in res if r["kind"] == kind and r["Lambda"] == lt]
            bp = max(rs, key=lambda r: r["F_plain"]); bi = max(rs, key=lambda r: r["F_improved"])
            print(f"Lambda={lt}: long mollifier best plain F={bp['F_plain']:+.4f} (y0={bp['y0']}) | best improved F={bi['F_improved']:+.4f} (y0={bi['y0']})")
    print(f"{time.time()-T:.0f}s")

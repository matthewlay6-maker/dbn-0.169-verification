"""Choose (y0, taper) for Lambda = 0.169 with the b-weighted long mollifier (y = 300, polynomial taper), float scouting.
Parallel over y0; each optimises (a2, a3) for the single-N bound F at N0, then reports the narrow-block bound."""
import json, math, sys, time
from multiprocessing import Pool
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, "levinson")
import levinson_cloud as L

LAM, Y = 0.169, 300


def run(y0):
    t0 = LAM - y0 * y0 / 2
    N0 = int(math.floor(math.sqrt(L.X / (4 * math.pi) + t0 / 16)))
    S = L.Lev(N0, t0, y0, Y)
    r = minimize(lambda x: -S.F(S.lam(x, "poly")), np.array([-1.7, 0.0]), method="Nelder-Mead",
                 options=dict(maxiter=120, xatol=1e-4, fatol=1e-7))
    lam = S.lam(r.x, "poly")
    hi, lo = S.Mbounds(lam)
    return dict(y0=y0, t0=t0, N0=N0, F=-r.fun, coeffs=list(r.x), M_hi=hi, M_lo=lo, evals=r.nfev)


if __name__ == "__main__":
    T = time.time()
    with Pool(5) as p:
        res = p.map(run, [0.06, 0.07, 0.08, 0.09, 0.10])
    for r in res:
        print(f"y0={r['y0']}: F={r['F']:+.5f} coeffs={[round(c, 4) for c in r['coeffs']]} max|M|={r['M_hi']:.3f} min|M|>={r['M_lo']:+.3f} ({r['evals']} evals)")
    json.dump(res, open("results/opt169.json", "w"), indent=1)
    print(f"{time.time()-T:.0f}s")

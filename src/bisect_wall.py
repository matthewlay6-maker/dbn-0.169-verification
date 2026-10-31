"""Bisect the lowest Lambda = t0 + y0^2/2 our method certifies-in-principle at RH height T_PT.
Known: 0.175 fails, 0.1787854 passes.  For each target: try y0 on a grid, t0 = Lambda - y0^2/2, at two
barrier locations (ours and Gomila's), and require, all with safety factors:
  middle  : first block [N0, 1.04 N0] (binding in our 0.18 run) memory-light F >= 10 * need
  tail    : A + B < 1.95 at some N1 <= 2e8
  barrier : min |f_t| on the boundary at t in {0, 0.05, 0.1, t0} >= 10 * e_barrier, winding 0
Proxy for the full run (the remaining blocks have so far been easier than the first)."""
import json, math, time
import numpy as np
import cost_map as CM
import barrier_gen as B
from gonogo_helpers import F_block_light
import full_run as FR

LOCS = {"ours": "6000347387366", "gomila": "6000000185827.5"}
YGRID = [0.12, 0.14, 0.16, 0.18, 0.20, 0.22]
OUT = "results/bisect_wall.json"


def tail_ok(t0, y0):
    for N1 in (2e7, 5e7, 1e8, 2e8):
        n = np.arange(1, 200001, dtype=float); ln = np.log(n)
        s1 = (1 + y0) / 2 + t0 / 2 * math.log(N1); s2 = s1 - y0
        A = np.sum(np.exp(t0/4*ln**2 - s1*ln)) + 1.006 * N1**-y0 * np.sum(np.exp(t0/4*ln**2 - s2*ln))
        mx = lambda s: max(2e5**(1-s)*math.exp(t0/4*math.log(2e5)**2), N1**(1-s)*math.exp(t0/4*math.log(N1)**2))
        if A + (mx(s1) + 1.006*N1**-y0*mx(s2)) * math.log(N1/2e5) < 1.95:
            return N1
    return None


def try_target(lam):
    best = None
    for loc, x0 in LOCS.items():
        X = float(x0) - 0.5
        for y0 in YGRID:
            t0 = lam - y0 * y0 / 2
            if t0 <= 0.05:
                continue
            N0 = int(math.floor(math.sqrt(X / (4 * math.pi) + t0 / 16)))
            need = CM.mollifier_max(t0, y0, (1 + y0) / 2 + t0 / 2 * math.log(N0)) * FR.e_bound(X, y0, t0)
            Fb = F_block_light(N0, int(N0 * 1.04) + 1, t0, y0)
            row = dict(loc=loc, y0=y0, t0=t0, F_block=Fb, need=need)
            if Fb < 10 * need:
                continue
            N1 = tail_ok(t0, y0)
            row["N1"] = N1
            if N1 is None:
                continue
            if best is None or Fb > best["F_block"]:
                best = row
    if best is None:
        return False, None
    # barrier samples for the best candidate
    B.configure(LOCS[best["loc"]], best["y0"], best["t0"])
    eb = FR.e_bound(float(LOCS[best["loc"]]) - 0.5, best["y0"], best["t0"], drop_t_term=True)
    ln, th = B.precompute(); pts, h0 = B.mesh(2 ** 12); pa = np.array(pts); G0, C = B.gamma_parts(pts)
    mins = []
    for ts in (0.0, 0.05, 0.1, best["t0"]):
        M1, M2, *_ = B.moments(ts, ln, th)
        f, _ = B.evaluate(pa, G0, C, ts, M1, M2)
        w = int(round(float(np.sum(np.angle(np.roll(f, -1) / f))) / (2 * math.pi)))
        mins.append((float(np.abs(f).min()), w))
    best.update(barrier_min=min(m for m, _ in mins), barrier_windings=[w for _, w in mins], e_barrier=eb)
    ok = best["barrier_min"] >= 10 * eb and all(w == 0 for _, w in mins)
    return ok, best


if __name__ == "__main__":
    lo, hi = 0.175, 0.1787854
    log = []
    t00 = time.time()
    for it in range(8):
        mid = round((lo + hi) / 2, 6)
        ok, best = try_target(mid)
        log.append(dict(Lambda=mid, ok=ok, best=best))
        print(f"[{it}] Lambda={mid:.6f}: {'PASS' if ok else 'fail'} {best if best else ''} [{time.time()-t00:.0f}s]", flush=True)
        if ok:
            hi = mid
        else:
            lo = mid
        json.dump(dict(lo_fail=lo, hi_pass=hi, log=log), open(OUT, "w"), indent=1)
    print(f"wall: fails at {lo:.6f}, passes at {hi:.6f}", flush=True)

"""Go/no-go for a certified Lambda bound below 0.2 using RH to Platt-Trudgian's height.
Hypothesis (i) of Theorem 1.2 needs X/2 <= T_PT = 3,000,175,332,800, so X <= 6,000,350,665,600.

1. Barrier location: X0 in the last 4e6 below the limit, maximising |prod_{p<=11} (1 - p^{-(1-iX)/2})^{-1}|
   (Polymath15's heuristic), with N constant on [X0-0.5, X0+0.5].
2. Middle range, hardest point: F(N0; t0, y0) (Euler-5 mollified, uniform in T, exact truncation at N0)
   against the needed |E|max * (e_A + e_B + e_C0); and the memory-light block bound on [N0, 1.002 N0]
   (fixed coefficients for m <= N-, worst case for N- < m <= N+).
3. Barrier at sample times: min |f_t| on the boundary of [X0 +- 0.5] x [y0, 1], derivative bounds,
   and an estimate of the number of time steps for the full run.
"""
import json, math, sys, time
import numpy as np
import flint
from flint import arb, acb
import cost_map as CM
import barrier_gen as B

XMAX = 6000350665600
OUT = "results/gonogo.json"
t00 = time.time()
res = {}

# ---- 1. barrier location
X = np.arange(XMAX - 4_000_000, XMAX - 1, dtype=np.float64) + 0.5   # candidate X0 (x = X0 +- 0.5 <= XMAX)
logabs = np.zeros_like(X)
for p in (2, 3, 5, 7, 11):
    ph = (X / 2 * math.log(p)) % (2 * math.pi)
    z = p ** -0.5 * np.exp(1j * ph)                 # p^{-(1-iX)/2} = p^{-1/2} e^{i (X/2) log p}
    logabs -= np.log(np.abs(1 - z))
order = np.argsort(-logabs)
cands = []
for j in order[:200]:
    x0 = int(X[j] - 0.5)
    try:
        N = B.configure(str(x0), 0.1, 0.2)
    except AssertionError:
        continue
    cands.append((x0, float(np.exp(logabs[j])), N))
    if len(cands) == 5:
        break
x0, eul, N0 = cands[0]
res["barrier_location"] = dict(X0=x0, euler_product=eul, N0=N0, candidates=cands)
print(f"1. barrier X0 = {x0} (Euler product {eul:.2f}), N0 = {N0}  [{time.time()-t00:.0f}s]", flush=True)
json.dump(res, open(OUT, "w"), indent=1)

# ---- 2. middle range at the hardest point
params = [(0.180, 0.1414), (0.175, 0.1732), (0.185, 0.1000), (0.175, 0.1414), (0.170, 0.1732)]
rows = []
for t0, y0 in params:
    Xf = float(x0)
    ok, F, need, N = CM.feasible(Xf, t0, y0)
    rows.append(dict(t0=t0, y0=y0, Lambda=t0 + y0 * y0 / 2, F_single_N0=F, need=need, ratio=F / need))
    print(f"2. t0={t0} y0={y0} Lambda={t0 + y0*y0/2:.4f}: F(N0)={F:.4f} need={need:.2e} {'OK' if ok else 'FAIL'}", flush=True)
res["middle_range"] = rows
json.dump(res, open(OUT, "w"), indent=1)


def F_block_light(Nm, Np, t0, y0):
    """Memory-light block bound: coefficients with m <= Nm fixed; Nm < m <= Np counted in absolute value."""
    divs = [1, 2, 3, 5, 6, 10, 15, 30]
    lam = {d: math.prod(-math.exp(t0 / 4 * math.log(p) ** 2) for p in (2, 3, 5) if d % p == 0) for d in divs}
    L = 30 * Np + 1
    bc = np.zeros(L); ac = np.zeros(L); var = np.zeros(L)
    m = np.arange(1, Np + 1, dtype=np.float64)
    lm = np.log(m)
    b = np.exp(t0 / 4 * lm * lm)
    a = 1.03 * np.exp(y0 * lm) * b
    core = m <= Nm
    c = Nm ** -y0
    for d in divs:
        idx = d * np.arange(1, Np + 1)
        bc[idx[core]] += lam[d] * b[core]
        ac[idx[core]] += lam[d] * a[core]
        var[idx[~core]] += abs(lam[d]) * (b[~core] + c * a[~core])
    sig = (1 + y0) / 2 - 0.001 + t0 / 2 * math.log(Nm)
    a1 = ac[1]
    r = (1 - c * a1) / (1 + c * a1)
    n = np.arange(L, dtype=np.float64)
    w = np.zeros(L); w[2:] = np.exp(-sig * np.log(n[2:]))
    num = np.maximum(np.abs(bc - c * ac), r * np.abs(bc + c * ac)) + var
    num[:2] = 0
    return 1 - c * a1 - float(np.sum(num * w))


best = max(rows, key=lambda r: (r["ratio"] > 1, -r["Lambda"]))
t0, y0 = best["t0"], best["y0"]
blk = []
for frac in (1.0005, 1.002):
    Fb = F_block_light(N0, int(N0 * frac), t0, y0)
    blk.append(dict(Nplus=int(N0 * frac), F_block=Fb))
    print(f"2b. block [N0, {frac}·N0] light bound at t0={t0}, y0={y0}: F = {Fb:.4f}", flush=True)
res["block_light"] = dict(t0=t0, y0=y0, rows=blk)
json.dump(res, open(OUT, "w"), indent=1)

# ---- 3. barrier samples
B.configure(str(x0), y0, t0)
ln, th = B.precompute()
pts, h0 = B.mesh(2 ** 14)
pts_arr = np.array(pts)
G0, C = B.gamma_parts(pts)
samp = []
for ts in (0.0, 0.03, 0.06, 0.1, 0.14, t0):
    Dz, Dt = B.deriv_bounds(ts, ln)
    M1, M2, a1, a2 = B.moments(ts, ln, th)
    f, _ = B.evaluate(pts_arr, G0, C, ts, M1, M2)
    mn = float(np.abs(f).min())
    wind = int(round(float(np.sum(np.angle(np.roll(f, -1) / f))) / (2 * math.pi)))
    need_h = 2 * 0.25 * mn / Dz
    samp.append(dict(t=ts, min_abs_f=mn, Dz=Dz, Dt=Dt, winding_on_mesh=wind, mesh_h=h0, h_needed=need_h,
                     est_dt=0.5 * mn / Dt))
    print(f"3. t={ts:.3f}: min|f|={mn:.3f} winding={wind} Dz={Dz:.3g} Dt={Dt:.3g} h_needed={need_h:.1e} "
          f"(mesh {h0:.1e}) dt~{0.5*mn/Dt:.1e}  [{time.time()-t00:.0f}s]", flush=True)
# step-count estimate: integrate dt over t
ts_ = np.array([s["t"] for s in samp]); dts = np.array([s["est_dt"] for s in samp])
grid = np.linspace(0, t0, 2001)
est_steps = float(np.trapezoid(1 / np.interp(grid, ts_, dts), grid))
res["barrier_samples"] = dict(rows=samp, est_steps=est_steps)
res["runtime_s"] = time.time() - t00
json.dump(res, open(OUT, "w"), indent=1)
print(f"estimated barrier time steps for t in [0,{t0}]: {est_steps:.0f}; go/no-go done in {time.time()-t00:.0f}s", flush=True)

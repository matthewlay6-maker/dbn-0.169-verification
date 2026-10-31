"""Go/no-go probe for lower targets at the same barrier location (results/gonogo.json)."""
import json, math, time, sys
import numpy as np
import cost_map as CM
import barrier_gen as B
from gonogo_helpers import F_block_light

g = json.load(open("results/gonogo.json"))
x0 = g["barrier_location"]["X0"]
cands = [(0.170, 0.1414), (0.165, 0.1732), (0.165, 0.1414), (0.160, 0.1732)]
out = []
t00 = time.time()
for t0, y0 in cands:
    lam = t0 + y0 * y0 / 2
    ok, F, need, N = CM.feasible(float(x0), t0, y0)
    Fb = F_block_light(N, int(N * 1.002), t0, y0)
    B.configure(str(x0), y0, t0)
    ln, th = B.precompute()
    pts, h0 = B.mesh(2 ** 13)
    pa = np.array(pts)
    G0, C = B.gamma_parts(pts)
    samp = []
    for ts in (0.0, 0.05, 0.1, t0):
        Dz, Dt = B.deriv_bounds(ts, ln)
        M1, M2, *_ = B.moments(ts, ln, th)
        f, _ = B.evaluate(pa, G0, C, ts, M1, M2)
        mn = float(np.abs(f).min())
        wind = int(round(float(np.sum(np.angle(np.roll(f, -1) / f))) / (2 * math.pi)))
        samp.append(dict(t=ts, min_abs_f=mn, winding=wind, Dt=Dt, est_dt=0.5 * mn / Dt))
    ts_ = np.array([s["t"] for s in samp]); dts = np.array([s["est_dt"] for s in samp])
    grid = np.linspace(0, t0, 2001)
    steps = float(np.trapezoid(1 / np.interp(grid, ts_, dts), grid))
    row = dict(t0=t0, y0=y0, Lambda=lam, F_N0=F, need=need, F_block=Fb,
               barrier_min_abs_f=min(s["min_abs_f"] for s in samp),
               barrier_windings=[s["winding"] for s in samp], est_steps=steps)
    out.append(row)
    print(f"Lambda={lam:.4f} (t0={t0}, y0={y0}): F(N0)={F:.4f} block={Fb:.4f} need={need:.1e} | "
          f"barrier min|f|={row['barrier_min_abs_f']:.3f} windings={row['barrier_windings']} steps~{steps:.0f} "
          f"[{time.time()-t00:.0f}s]", flush=True)
json.dump(out, open("results/gonogo2.json", "w"), indent=1)

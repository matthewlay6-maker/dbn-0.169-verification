"""Session 3A: spot-check of the exterior claim in the 0.1595 certificate.
Claim (as stated in its ASSEMBLY.md, read-only): for t in [319/2000, 3/16] and all x >= X = 6000000185827,
(H_t(x), H_t'(x)) != (0,0), via |u_x| > 1/400000 at real zeros, u = H_t/|M_t|.

At a real zero, |u_x| = |d/dx (H_t/B_t)| (since |B_t| = |M_t(s_+)|), and H_t/B_t = f_t + O(e), e from Theorem 1.3.
We sample windows of x near X and at larger heights, locate sign changes of Re(f_t * phase) (H_t is real on the
real axis, so H_t/B_t = |H_t/B_t| e^{-i arg B_t}), and record |f_t'| at each zero (finite difference).
This is a plausibility probe only: a sample can find a counterexample, never confirm 'for all x'.
float64 with phases reduced mod 2 pi in arb; errors are reported, not certified.
"""
import json, math, sys, time
import numpy as np
import flint
from flint import arb, acb
import pm15

flint.ctx.prec = 128
PI = arb.pi()
OUT = "results/audit_spotcheck.json"


def window(Xc, length, t, step=0.02):
    """f_t on [Xc, Xc+length], y = 0, using Taylor-free direct sums with reduced phases."""
    Xa = arb(repr(Xc))
    N = int(((Xa / (4 * PI) + arb(repr(t)) / 16).sqrt()).floor().unique_fmpz())
    n = np.arange(1, N + 1, dtype=np.float64)
    ln = np.log(n)
    th = np.empty(N)
    half = Xa / 2
    tp = 2 * PI
    for k in range(1, N + 1):
        v = half * arb(k).log()
        th[k - 1] = float(v - tp * (v / tp).floor())
    # reference quantities at Xc (vary by O(length/X) across the window: allowance below)
    P = pm15.parts(arb(repr(t)), Xa, arb(0))
    a1 = complex(P["sstar"] - (1 - acb(0, 1) * Xa) / 2)   # (t/2) alpha
    b = np.exp(t / 4 * ln * ln)
    base1 = b * np.exp(-0.5 * ln - a1 * ln) * np.exp(1j * th)
    # second sum: n^{-s**}, s** = (1+ix)/2 + (t/2) alpha((1+ix)/2) = conj(s*) (y = 0)
    base2 = b * np.exp(-0.5 * ln - np.conj(a1) * ln) * np.exp(-1j * th)
    xs = np.arange(0, length, step)
    vals = np.empty(len(xs), complex)
    tt = arb(repr(t))
    for j, dx in enumerate(xs):
        x = Xa + arb(repr(float(dx)))
        gam = complex(pm15.Mt(tt, (1 + acb(0, 1) * x) / 2) / pm15.Mt(tt, (1 - acb(0, 1) * x) / 2))
        ph = np.exp(0.5j * dx * ln)
        vals[j] = np.sum(base1 * ph) + gam * np.sum(base2 / ph)
    return xs, vals, N


def zeros_and_slopes(xs, vals):
    # zeros of H_t on the real axis = zeros of |H_t/B_t| ~ |f_t|: take local minima of |f| that are
    # consistent with a crossing inside one grid step (|f| < 0.6 * slope * step)
    out = []
    mag = np.abs(vals)
    for j in range(1, len(xs) - 1):
        if mag[j] <= mag[j - 1] and mag[j] <= mag[j + 1]:
            h = xs[1] - xs[0]
            slope = abs(vals[j + 1] - vals[j - 1]) / (2 * h)
            if mag[j] < 0.6 * slope * h:
                out.append((float(xs[j]), float(mag[j]), float(slope)))
    return out


def main():
    t = 319 / 2000
    X = 6000000185827.0
    centres = [X, X + 1e6, X + 1e8, 2e13, 4e13, 6.6e13]
    length = float(sys.argv[1]) if len(sys.argv) > 1 else 200.0
    res = dict(t=t, windows=[])
    t00 = time.time()
    for Xc in centres:
        t0 = time.time()
        xs, vals, N = window(Xc, length, t)
        mins = zeros_and_slopes(xs, vals)
        lx = math.log(Xc / (4 * math.pi))
        ec0 = (Xc / (4 * math.pi)) ** (-0.25) * math.exp(-t / 16 * lx * lx + 1.24 * 2 / (N - 0.125)
                                                         + (3 * abs(complex(lx, math.pi / 2)) + 10.44) / (Xc - 12))
        near_zero = mins
        res["windows"].append(dict(centre=Xc, length=length, N=N, n_local_minima=len(mins),
                                   n_near_zeros=len(near_zero),
                                   min_slope_at_near_zero=min((m[2] for m in near_zero), default=None),
                                   median_slope=float(np.median([m[2] for m in near_zero])) if near_zero else None,
                                   eC0_bound=ec0, runtime_s=time.time() - t0))
        w = res["windows"][-1]
        print(f"x~{Xc:.4g}: N={N} zeros~{w['n_near_zeros']} min|f'| at zeros={w['min_slope_at_near_zero']} "
              f"median={w['median_slope']} eC0<={ec0:.2e} {time.time()-t0:.0f}s", flush=True)
        json.dump(res, open(OUT, "w"), indent=1)
    res["runtime_s"] = time.time() - t00
    json.dump(res, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main()

"""Session 2b, claim (a): Polymath15 barrier at X0 = 6e10 + 83952.
Region: X0-0.5 <= x <= X0+0.5, 0.2 <= y <= 1, 0 <= t <= 0.2, N = N0 = 69098 (constant).
Show: for every t in [0, 0.2], f_t(x+iy) stays outside |z| <= 1.25e-3 on the rectangle boundary
and winds zero times around 0 (then by Rouche + the sweep proposition, H_t has no zero in R).

Method (independent implementation of the paper's algorithm):
  f_t = S1 + gamma*S2,  S1 = sum b_n n^{-s*},  S2 = sum b_n n^{-s**}  (paper eq. 'ftne'),
  s* = (1+y-ix)/2 + (t/2) alpha((1+y-ix)/2),  s** = (1-y+ix)/2 + (t/2) alpha((1-y+ix)/2).
  Around the centre (x, y) = (X0, 0.6): n^{-s*} = n^{-s*_c} n^{w},  n^{-s**} = n^{-s**_c} n^{-w},
  w = (i(x-X0) - (y-0.6))/2, |w| <= 0.32, up to the alpha-variation allowance.
  S1, S2 by a 40-term Taylor expansion in w (moments sum a_n log^k n), with explicit tail bound.
  Large phases (X0/2) log n are reduced mod 2 pi in arb.  gamma = G0 * exp(t c), with G0, c from arb.
  Mesh on the boundary: spacing h; derivative bounds D_z, D_t from the paper's Lemma (valid for all
  t' >= t, since b_n^t n^{-(t/4) log(x/4pi)} is non-increasing in t).  Step:
     |f_{t'}(z) - f_t(z_j)| <= D_z h/2 + D_t (t'-t) + allowances  =: rho
  and we require min_j |f_t(z_j)| - rho > 1.25e-3; the winding number of the polygon (arg increments
  each < pi) then equals that of the true curve for all t' in [t, t+dt].
  NOTE: the paper writes D_t|t'-t|/(2n); we use D_t|t'-t| (the t-derivative is not reduced by mesh density).
"""
import json, math, time, sys
import numpy as np
import flint
from flint import arb, acb

flint.ctx.prec = 128
X0s = "60000083952"
X0 = arb(X0s)
N = 69098
Y0, Y1, YC = 0.2, 1.0, 0.6
TMAX = 0.2


def configure(x0s, y0, tmax):
    """Re-target the module: barrier centre x0s (string), canopy height y0, final time tmax."""
    global X0s, X0, N, Y0, YC, TMAX
    X0s, X0 = x0s, arb(x0s)
    Y0, YC, TMAX = y0, (y0 + 1.0) / 2, tmax
    lo = int(((X0 - arb("0.5")) / (4 * PI)).sqrt().floor().unique_fmpz())
    hi = int((((X0 + arb("0.5")) / (4 * PI) + arb(repr(tmax)) / 16).sqrt()).floor().unique_fmpz())
    assert lo == hi, "N not constant on the barrier"
    N = lo
    return N
K = 40
MARGIN = 1.25e-3
OUT = "results/barrier_pm15.json"
PI = arb.pi()


def alpha(s):
    return 1 / (2 * s) + 1 / (s - 1) + (s / (2 * PI)).log() / 2


def logM0(s):
    return s.log() + (s - 1).log() - s / 2 * PI.log() + ((2 * PI).sqrt() / 16).log() \
        + (s / 2 - acb("0.5")) * (s / 2).log() - s / 2


def mesh(M):
    """Boundary points (delta, y), counter-clockwise, spacing <= 1/M."""
    pts = []
    nb = M
    hgt = Y1 - Y0
    ny = int(math.ceil(hgt * M))
    for k in range(nb):          # bottom: y = Y0, delta -0.5 -> 0.5
        pts.append((-0.5 + k / nb, Y0))
    for k in range(ny):          # right: delta = 0.5, y Y0 -> 1
        pts.append((0.5, Y0 + hgt * k / ny))
    for k in range(nb):          # top: y = 1, delta 0.5 -> -0.5
        pts.append((0.5 - k / nb, Y1))
    for k in range(ny):          # left: delta = -0.5, y 1 -> Y0
        pts.append((-0.5, Y1 - hgt * k / ny))
    return pts, 1.0 / M


def gamma_parts(pts):
    """For each mesh point: G0 = M0(s)/M0(1-s) (as complex) and c = (alpha(s)^2 - alpha(1-s)^2)/4,
    s = (1-y+ix)/2, so gamma_t = G0 * exp(t c)."""
    G0 = np.empty(len(pts), complex)
    C = np.empty(len(pts), complex)
    for j, (d, y) in enumerate(pts):
        x = X0 + arb(repr(d))
        s = (1 - arb(repr(y)) + acb(0, 1) * x) / 2
        g = (logM0(s) - logM0(1 - s)).exp()
        c = (alpha(s) ** 2 - alpha(1 - s) ** 2) / 4
        G0[j] = complex(g)
        C[j] = complex(c)
    return G0, C


def precompute():
    n = np.arange(1, N + 1, dtype=np.float64)
    ln = np.log(n)
    th = np.empty(N)
    half = X0 / 2
    tp = 2 * PI
    for k in range(1, N + 1):
        v = half * arb(k).log()
        v = v - tp * (v / tp).floor()
        th[k - 1] = float(v)
    return ln, th


def moments(t, ln, th):
    """Moments for S1 and S2 at the centre, for time t."""
    xc = acb(arb(1) + arb(repr(YC)), -X0) / 2       # (1+y-ix)/2 at centre
    xcc = acb(arb(1) - arb(repr(YC)), X0) / 2       # (1-y+ix)/2 at centre
    a1 = complex(arb(repr(t)) / 2 * alpha(xc))
    a2 = complex(arb(repr(t)) / 2 * alpha(xcc))
    b = np.exp(t / 4 * ln * ln)
    # n^{-s*_c} = exp(-(1+YC)/2 ln) * exp(+i (X0/2) ln) * exp(-a1 ln)
    c1 = b * np.exp(-(1 + YC) / 2 * ln - a1 * ln) * np.exp(1j * th)
    c2 = b * np.exp(-(1 - YC) / 2 * ln - a2 * ln) * np.exp(-1j * th)
    M1 = np.empty(K + 1, complex)
    M2 = np.empty(K + 1, complex)
    p = np.ones_like(ln)
    for k in range(K + 1):
        M1[k] = np.sum(c1 * p)
        M2[k] = np.sum(c2 * p)
        p = p * ln
    abs1 = float(np.sum(np.abs(c1) * np.exp(0.34 * ln)))
    abs2 = float(np.sum(np.abs(c2) * np.exp(0.34 * ln)))
    return M1, M2, abs1, abs2


def deriv_bounds(t, ln):
    """D_z, D_t upper bounds (paper's Lemma), valid for all t' in [t, 0.2] on the boundary."""
    x = float(X0s) - 0.5
    lx = math.log(x / (4 * math.pi))
    sig = (1 + Y0) / 2 + t / 4 * lx - 1e-18
    w = np.exp(t / 4 * ln * ln - sig * ln)          # b_n / n^{Re s*} (worst case y = 0.2)
    g = math.exp(0.02) * 1.000001                   # |gamma| n^y N^{|kappa|} <= e^{0.02} (n^2 <= x/4pi)
    assert N * N <= x / (4 * math.pi)
    big = math.log(math.hypot(2.0, x) / (4 * math.pi)) + math.pi + 3 / x
    Dz = np.sum(w * (ln / 2 + t * ln / (4 * (x - 6)))) \
        + g * np.sum(w * (t * ln / (4 * (x - 6)) + big * (0.5 + t / (4 * (x - 6)))))
    q = 0.25 * ln * np.log(x / (4 * math.pi) / np.exp(ln)) + math.pi / 8 * ln + 2 * ln / (x - 6)
    Dt = np.sum(w * q) + g * np.sum(w * (q + 0.25 * (math.pi / 2 + 8 / (x - 6)) * (lx + 8 / (x - 6))))
    return float(Dz) * 1.0001, float(Dt) * 1.0001


def evaluate(pts_arr, G0, C, t, M1, M2):
    d = pts_arr[:, 0]
    y = pts_arr[:, 1]
    w = (1j * d - (y - YC)) / 2
    S1 = np.zeros(len(d), complex)
    S2 = np.zeros(len(d), complex)
    wk = np.ones(len(d), complex)
    fact = 1.0
    for k in range(K + 1):
        S1 += wk * M1[k] / fact
        S2 += ((-1) ** k) * wk * M2[k] / fact
        wk = wk * w
        fact *= (k + 1)
    gam = G0 * np.exp(t * C)
    return S1 + gam * S2, np.abs(gam).max()


def main():
    t_start = time.time()
    ln, th = precompute()
    lev = 15
    pts, h0 = mesh(2 ** lev)
    pts_arr = np.array(pts)
    G0, C = gamma_parts(pts)
    print(f"precompute done ({len(pts)} fine mesh points) {time.time()-t_start:.0f}s", flush=True)
    L = math.log(N)
    # Taylor tail (|w| <= 0.34 for y0 >= 0.1): (0.34 L)^{K+1}/(K+1)! e^{0.33 L} per unit |a_n|
    tail_rel = (0.34 * L) ** (K + 1) / math.factorial(K + 1) * math.exp(0.34 * L)
    t = 0.0
    steps = []
    while t < TMAX:
        Dz, Dt = deriv_bounds(t, ln)
        M1, M2, abs1, abs2 = moments(t, ln, th)
        # allowances: Taylor tail, float rounding in moments (N eps), alpha/kappa variation (2e-11 rel)
        allow = (tail_rel + 2e-11 + N * 2.3e-16 * 4) * (abs1 + math.exp(0.02) * abs2)
        # choose mesh level: largest spacing with D_z h/2 <= 0.2 * (estimated min |f|)
        f_coarse, _ = evaluate(pts_arr[:: 2 ** 6], G0[:: 2 ** 6], C[:: 2 ** 6], t, M1, M2)
        mest = np.abs(f_coarse).min()
        step = 2 ** 6
        while step > 1 and Dz * (h0 * step) / 2 > 0.2 * mest:
            step //= 2
        P = pts_arr[::step]
        f, gmax = evaluate(P, G0[::step], C[::step], t, M1, M2)
        mn = float(np.abs(f).min())
        h = h0 * step
        spare = mn - MARGIN - Dz * h / 2 - allow
        if spare <= 0:
            print(f"FAIL at t={t}: min|f|={mn:.4e}, Dz h/2={Dz*h/2:.3e}", flush=True)
            steps.append(dict(t=t, fail=True, min_abs_f=mn, Dz=Dz, h=h))
            break
        ang = np.angle(np.roll(f, -1) / f)
        assert np.max(np.abs(ang)) < math.pi * 0.99
        wind = int(round(float(np.sum(ang)) / (2 * math.pi)))
        dt = spare / Dt
        steps.append(dict(t=t, dt=dt, points=len(P), h=h, min_abs_f=mn, Dz=Dz, Dt=Dt,
                          allowance=allow, winding=wind))
        if wind != 0:
            print(f"winding {wind} at t={t}", flush=True)
            break
        if len(steps) % 50 == 1:
            print(f"t={t:.6f} pts={len(P)} min|f|={mn:.3f} Dz={Dz:.3g} Dt={Dt:.3g} dt={dt:.2e} "
                  f"{time.time()-t_start:.0f}s", flush=True)
            json.dump(dict(steps=steps), open(OUT, "w"))
        t = t + dt
    ok = all(not s.get("fail") and s.get("winding") == 0 for s in steps) and t >= TMAX
    res = dict(ok=ok, n_steps=len(steps), t_reached=t, min_abs_f_overall=min(s["min_abs_f"] for s in steps),
               max_points=max(s.get("points", 0) for s in steps), runtime_s=time.time() - t_start,
               steps=steps)
    json.dump(res, open(OUT, "w"))
    print(f"barrier (a): ok={ok} steps={len(steps)} min|f|={res['min_abs_f_overall']:.4f} "
          f"max mesh={res['max_points']} {time.time()-t_start:.0f}s", flush=True)


if __name__ == "__main__":
    main()

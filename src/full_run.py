"""Full computation toward Lambda <= t0 + y0^2/2 = 0.18 (t0 = 0.17, y0 = 0.1414) via Polymath15's Theorem 1.2,
with RH verified to T_PT = 3,000,175,332,800 (Platt-Trudgian 2021) as the only external input.

Structure follows Polymath15 Section 8 (our reproduction in repro_pm15.py / barrier_pm15.py):
  (i)   RH to X/2: X = X0 - 0.5 <= 2 T_PT                                  [external input]
  (iii) barrier [X0-0.5, X0+0.5] x [y0, 1], t in [0, t0]: winding 0, |f| > e  [barrier_gen]
  (ii)  region x >= X, y0 <= y <= 1, t = t0, by the argument principle for E_{t,5} H/B on
        R = {N0 <= N <= N1}: lower edge (mollified block bounds), upper edge y = 1, right edge N = N1
        (tail bound), left edge x = X0 - 0.5 (barrier output at t = t0); and N >= N1 by the tail bound.
Arithmetic: arb for constants and error bounds; float64 with explicit allowances for the long sums.
Usage: python full_run.py [errors|barrier|middle|tail|edges|all]
"""
import json, math, sys, time, os
import numpy as np
import flint
from flint import arb, acb

T0, Y0 = 0.17, 0.1414
TPT = 3000175332800
OUT = os.environ.get("OUT", "results/full_run.json")
X0 = json.load(open("results/gonogo.json"))["barrier_location"]["X0"]
X = X0 - 0.5
assert X <= 2 * TPT
N1 = 50_000_000
PI = arb.pi()
flint.ctx.prec = 128
N0 = int(((arb(repr(X)) / (4 * PI) + arb(repr(T0)) / 16).sqrt()).floor().unique_fmpz())


def load():
    try:
        return json.load(open(OUT))
    except FileNotFoundError:
        return {}


def save(d):
    json.dump(d, open(OUT, "w"), indent=1)


def e_bound(x, y, t, drop_t_term=False):
    """Upper bound for e_A + e_B + e_{C,0} at (x, y, t) (Theorem 1.3 / Prop. 'estimates'), arb.
    Monotone decreasing in x and y (Polymath15 Lemma elem-lem(vi)); with drop_t_term the
    -t/16 log^2 term is discarded so the bound holds for all 0 <= t' <= t."""
    x, y, t = arb(repr(x)), arb(repr(y)), arb(repr(t))
    lx = (x / (4 * PI)).log()
    N = ((x / (4 * PI) + t / 16).sqrt()).floor()
    d1 = (t * t / 16 * lx * lx + arb("0.626")) / (x - arb("6.66"))
    # e_A + e_B <= (e^{d1}-1)(1 + 1.021) F, F <= N^{(1-y)/2 ... } use sum_{n<=N} n^{-(1+y)/2} <= N^{(1-y)/2}/((1-y)/2)
    s = (1 - y) / 2
    eab = (d1.exp() - 1) * (1 + arb("0.02").exp() * arb("1.0000001")) * (N + 1) ** s / s
    tt = arb(0) if drop_t_term else t
    ec0 = (-(1 + y) / 4 * lx).exp() * (-tt / 16 * lx * lx + (3 * abs(acb(lx, PI / 2)) + arb("3.58")) / (x - arb("8.52"))).exp() \
        * (1 + arb("1.24") * (arb(3) ** y + arb(3) ** (-y)) / (N - arb("0.125")) + arb("6.92") / (x - arb("6.66")))
    return float((eab + ec0).upper())


def mollifier(t, sig):
    """|E_{t,5}| upper/lower bounds and max |arg| for Re s >= sig."""
    hi, lo, ar = 1.0, 1.0, 0.0
    for p in (2, 3, 5):
        v = math.exp(t / 4 * math.log(p) ** 2) / p ** sig
        hi *= 1 + v; lo *= 1 - v; ar += math.asin(v)
    return hi * (1 + 1e-12), lo * (1 - 1e-12), ar * (1 + 1e-12)


def part_errors(d):
    eb_barrier = e_bound(X, Y0, T0, drop_t_term=True)       # all t in [0, t0], y >= y0, x >= X
    eb_t0 = e_bound(X, Y0, T0)                              # t = t0 region (ii)
    sig_min = (1 + Y0) / 2 + T0 / 2 * math.log(N0) - 1e-6   # Re s_* lower bound on region (ii)
    Ehi, Elo, Earg = mollifier(T0, sig_min)
    d["errors"] = dict(X0=X0, X=X, N0=N0, N1=N1, t0=T0, y0=Y0, Lambda=T0 + Y0 * Y0 / 2,
                       e_barrier=eb_barrier, e_t0=eb_t0, sig_min=sig_min,
                       E_hi=Ehi, E_lo=Elo, E_arg=Earg, need_mollified=Ehi * eb_t0)
    print("errors:", d["errors"], flush=True)


def part_barrier(d):
    import barrier_gen as B
    B.configure(str(X0), Y0, T0)
    B.MARGIN = d["errors"]["e_barrier"]
    B.OUT = "results/full_barrier.json"
    prev = json.load(open(B.OUT)) if os.path.exists(B.OUT) else {}
    if not (prev.get("ok") and abs(prev.get("t_reached", 0)) >= T0 and prev.get("n_steps")):
        B.main()
    b = json.load(open(B.OUT))
    d["barrier"] = {k: v for k, v in b.items() if k != "steps"}
    # left edge of region (ii) = barrier left side at t = t0: record |f| and arg there
    ln, th = B.precompute()
    pts, h0 = B.mesh(2 ** 13)
    left = np.array([p for p in pts if p[0] == -0.5])
    G0, C = B.gamma_parts([(float(p[0]), float(p[1])) for p in left])
    M1, M2, *_ = B.moments(T0, ln, th)
    f, _ = B.evaluate(left, G0, C, T0, M1, M2)
    Dz, _ = B.deriv_bounds(T0, ln)
    hy = (1 - Y0) / (len(left))
    d["left_edge"] = dict(min_abs_f=float(np.abs(f).min()), max_abs_arg=float(np.abs(np.angle(f)).max()),
                          Dz=Dz, mesh_h=hy, slack=Dz * hy / 2)
    print("barrier:", d["barrier"]["ok"], "left edge:", d["left_edge"], flush=True)


# ---- middle range: chunked, memory-light block bounds (lower edge y = y0) ----
DIVS = [1, 2, 3, 5, 6, 10, 15, 30]


def block_F(args):
    Nm, Np, t0, y0 = args
    lam = {d: math.prod(-math.exp(t0 / 4 * math.log(p) ** 2) for p in (2, 3, 5) if d % p == 0) for d in DIVS}
    c = Nm ** -y0
    sig = (1 + y0) / 2 - 0.001 + t0 / 2 * math.log(Nm)
    a1 = 1.03                     # alpha_1 = 1.03 * 1^{y0} * b_1
    r = (1 - c * a1) / (1 + c * a1)
    L = 30 * Np
    tot = 0.0
    CH = 4_000_000
    for s in range(2, L + 1, CH):
        e = min(L + 1, s + CH)
        n = np.arange(s, e, dtype=np.int64)
        bc = np.zeros(e - s); ac = np.zeros(e - s); var = np.zeros(e - s)
        for d in DIVS:
            mlo = -(-s // d); mhi = (e - 1) // d
            if mhi < mlo:
                continue
            m = np.arange(mlo, min(mhi, Np) + 1, dtype=np.float64)
            if len(m) == 0:
                continue
            lm = np.log(m)
            b = np.exp(t0 / 4 * lm * lm)
            a = 1.03 * np.exp(y0 * lm) * b
            pos = (m * d).astype(np.int64) - s
            core = m <= Nm
            # positions m*d - s are distinct within one divisor pass, so plain indexed adds are exact
            bc[pos[core]] += lam[d] * b[core]
            ac[pos[core]] += lam[d] * a[core]
            var[pos[~core]] += abs(lam[d]) * (b[~core] + c * a[~core])
        w = np.exp(-sig * np.log(n.astype(np.float64)))
        num = np.maximum(np.abs(bc - c * ac), r * np.abs(bc + c * ac)) + var
        tot += float(np.sum(num * w))
    F = 1 - c * a1 - tot
    # kappa term Z (paper: n^{|kappa|} - 1 with |kappa| <= t y/(2(x_N - 6))): bound crudely
    Z = 1.644 * 1e-9
    return dict(Nminus=Nm, Nplus=Np, F=F, F_lower=F - 1e-6 - Z, sigma=sig)


def part_middle(d, ratio=1.04, workers=6):
    from multiprocessing import Pool
    need = d["errors"]["need_mollified"]
    edges = [N0]
    while edges[-1] < N1:
        edges.append(min(N1, int(edges[-1] * ratio) + 1))
    jobs = [(a, b, T0, Y0) for a, b in zip(edges[:-1], edges[1:])]
    done = {(r["Nminus"], r["Nplus"]): r for r in d.get("middle", {}).get("rows", [])}
    todo = [j for j in jobs if (j[0], j[1]) not in done]
    print(f"middle: {len(jobs)} blocks, {len(todo)} to do, need F >= {need:.2e}", flush=True)
    t0 = time.time()
    with Pool(workers) as pool:
        for r in pool.imap_unordered(block_F, sorted(todo, key=lambda j: -j[1])):
            done[(r["Nminus"], r["Nplus"])] = r
            rows = sorted(done.values(), key=lambda r: r["Nminus"])
            d["middle"] = dict(rows=rows, need=need, n_blocks=len(jobs),
                               min_F=min(x["F_lower"] for x in rows),
                               ok=len(rows) == len(jobs) and all(x["F_lower"] >= need for x in rows))
            save(d)
            print(f"  block [{r['Nminus']},{r['Nplus']}] F={r['F']:.4f}  ({len(rows)}/{len(jobs)}, {time.time()-t0:.0f}s)", flush=True)
    print("middle ok:", d["middle"]["ok"], "min F:", d["middle"]["min_F"], flush=True)


def part_tail(d):
    """N >= N1, t = t0, y >= y0: |f_t - 1| <= A + B - 1 (Polymath15 Lemma 'largen', twice), arb for A."""
    t, y = arb(repr(T0)), arb(repr(Y0))
    N0l = 200_000
    lN = arb(N1).log()
    xN1 = 4 * PI * arb(N1) ** 2 - PI * t / 4
    corr = arb("1e-20")
    s1 = (1 + y) / 2 + t / 2 * lN - corr
    s2 = s1 - y - arb("1e-15")
    A1 = arb(0); A2 = arb(0)
    for n in range(1, N0l + 1):
        ln = arb(n).log()
        b = (t / 4 * ln * ln).exp()
        A1 += b * (-s1 * ln).exp()
        A2 += b * (-s2 * ln).exp()
    A = A1 + arb("1.006") * (xN1 / (4 * PI)) ** (-y / 2) * A2

    def Bof(N):
        N = arb(N); lN = N.log()
        s1 = (1 + y) / 2 + t / 2 * lN - corr
        s2 = s1 - y - arb("1e-15")
        bN0 = (t / 4 * arb(N0l).log() ** 2).exp(); bN = (t / 4 * lN ** 2).exp()
        m1 = max(arb(N0l) ** (1 - s1) * bN0, N ** (1 - s1) * bN, key=float)
        m2 = max(arb(N0l) ** (1 - s2) * bN0, N ** (1 - s2) * bN, key=float)
        return (m1 + arb("1.006") * N ** (-y) * m2) * (N / N0l).log()
    Bs = [float(Bof(int(N1 * 10 ** (k / 10))).upper()) for k in range(0, 351)]
    S = float(A.upper()) + max(Bs) - 1          # |f_t - 1| <= S for all N >= N1 (A decreasing in N)
    e = e_bound(4 * math.pi * N1 ** 2, Y0, T0)
    Ehi, Elo, Earg = d["errors"]["E_hi"], d["errors"]["E_lo"], d["errors"]["E_arg"]
    phi = math.asin(min(S, 1.0)) + Earg
    r = Elo * (1 - S)
    dist = r if phi <= math.pi / 2 else r * math.sin(math.pi - phi)
    d["tail"] = dict(N1=N1, A=float(A.upper()), B_max=max(Bs), B_decreasing=all(Bs[i + 1] <= Bs[i] for i in range(len(Bs) - 1)),
                     abs_f_minus_1_le=S, f_lower=1 - S, e=e, ok_region=bool(1 - S > e),
                     right_edge_dist=dist, right_edge_need=Ehi * e, ok_right_edge=bool(dist > Ehi * e and phi < math.pi))
    print("tail:", d["tail"], flush=True)


def part_edges(d):
    """Upper edge y = 1 (N0 <= N <= N1) and assembly of the argument-principle conditions."""
    t = T0
    sig1 = 1.0 + t / 2 * math.log(N0) - 1e-6            # Re s_* >= (1+1)/2 + (t/4) log(x/4pi)
    n = np.arange(2, N1 + 1, dtype=np.float64)
    ln = np.log(n)
    b = np.exp(t / 4 * ln * ln)
    s_first = float(np.sum(b * np.exp(-sig1 * ln)))      # sum_{n>=2} b_n n^{-sigma}
    # gamma term at y = 1: |gamma| n^y <= 1.03 n / N (|gamma| <= e^{0.02} (x/4pi)^{-1/2}, N <= sqrt(x/4pi + t/16)).
    # For n <= N0 use n/N <= n/N0; for N0 < n <= N use n/N <= 1.  Sum of both pieces bounds every N in [N0, N1].
    nn = np.arange(1, N1 + 1, dtype=np.float64); l2 = np.log(nn)
    bb = np.exp(t / 4 * l2 ** 2)
    s_second = 1.03 / N0 * float(np.sum(bb[:N0] * np.exp((1 - sig1) * l2[:N0])))         + 1.03 * float(np.sum(bb[N0:] * np.exp(-sig1 * l2[N0:])))
    S = (s_first + s_second) * (1 + 1e-9)
    Ehi, Elo, Earg = d["errors"]["E_hi"], d["errors"]["E_lo"], d["errors"]["E_arg"]
    phi = math.asin(min(S, 1)) + Earg
    r = Elo * (1 - S)
    dist = r if phi <= math.pi / 2 else r * math.sin(math.pi - phi)
    need = d["errors"]["need_mollified"]
    d["upper_edge"] = dict(abs_f_minus_1_le=S, dist=dist, need=need, ok=bool(S < 1 and dist > need and phi < math.pi))
    le = d.get("left_edge")
    if le:
        fl = le["min_abs_f"] - le["slack"]
        ar = le["max_abs_arg"] + math.asin(min(1, le["slack"] / le["min_abs_f"]))
        phi = ar + Earg
        r = Elo * fl
        dist = r if phi <= math.pi / 2 else r * math.sin(math.pi - phi)
        d["left_edge"].update(dist=dist, need=need, ok=bool(dist > need and phi < math.pi))
    checks = dict(barrier=d.get("barrier", {}).get("ok"), middle=d.get("middle", {}).get("ok"),
                  tail_region=d.get("tail", {}).get("ok_region"), right_edge=d.get("tail", {}).get("ok_right_edge"),
                  upper_edge=d["upper_edge"]["ok"], left_edge=d.get("left_edge", {}).get("ok"))
    d["verdict"] = dict(checks=checks, all_ok=all(v is True for v in checks.values()),
                        Lambda=T0 + Y0 * Y0 / 2, external_input="RH verified to 3,000,175,332,800 (Platt-Trudgian 2021)")
    print("edges/verdict:", d["upper_edge"], d.get("left_edge"), d["verdict"], flush=True)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    d = load()
    for name, fn in (("errors", part_errors), ("barrier", part_barrier), ("middle", part_middle),
                     ("tail", part_tail), ("edges", part_edges)):
        if which in (name, "all"):
            t0 = time.time()
            fn(d)
            d.setdefault("runtimes", {})[name] = time.time() - t0
            save(d)

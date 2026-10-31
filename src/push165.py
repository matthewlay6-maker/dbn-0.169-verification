"""Scouting (float64): how far can sharper constants and bigger Euler mollifiers push the lower-edge wall?
Single-N bound F(N0) at the hardest point, uniform in T (improved triangle inequality of Polymath15 section 8),
with mollifier E_{t,P} = prod_{p in P} (1 - b_p p^{-s}) and parameters:
  sig_slack : sigma = (1+y)/2 + (t/2) log N - sig_slack      (Polymath: 0.001; rigorous minimum ~1e-12)
  gcon      : |gamma| n^y <= gcon (n/N)^y                       (Polymath: 1.03; rigorous e^{0.02 y}(1+1e-9))
"""
import itertools, json, math, sys, time
import numpy as np

X = 6000347387366 - 0.5


def divisors(P):
    out = []
    for k in range(len(P) + 1):
        for c in itertools.combinations(P, k):
            out.append(math.prod(c))
    return sorted(out)


def F_single(N, t0, y0, P, sig_slack, gcon):
    D = math.prod(P)
    divs = divisors(P)
    lam = {d: math.prod(-math.exp(t0 / 4 * math.log(p) ** 2) for p in P if d % p == 0) for d in divs}
    sig = (1 + y0) / 2 + t0 / 2 * math.log(N) - sig_slack
    c = N ** -y0
    a1 = gcon
    r = (1 - c * a1) / (1 + c * a1)
    m = np.arange(1, N + 1, dtype=np.float64)
    lm = np.log(m)
    b = np.exp(t0 / 4 * lm * lm)
    a = gcon * np.exp(y0 * lm) * b
    L = D * N
    tot = 0.0
    CH = 20_000_000
    for s in range(2, L + 1, CH):
        e = min(L + 1, s + CH)
        bc = np.zeros(e - s); ac = np.zeros(e - s)
        for d in divs:
            mlo = -(-s // d); mhi = min((e - 1) // d, N)
            if mhi < mlo:
                continue
            idx = np.arange(mlo, mhi + 1) * d - s
            bc[idx] += lam[d] * b[mlo - 1:mhi]
            ac[idx] += lam[d] * a[mlo - 1:mhi]
        n = np.arange(s, e, dtype=np.float64)
        w = np.exp(-sig * np.log(n))
        tot += float(np.sum(np.maximum(np.abs(bc - c * ac), r * np.abs(bc + c * ac)) * w))
    return 1 - c * a1 - tot


if __name__ == "__main__":
    N0 = int(math.floor(math.sqrt(X / (4 * math.pi) + 0.17 / 16)))
    configs = [("Euler-5, Polymath constants", (2, 3, 5), 0.001, None),
               ("Euler-5, tight constants", (2, 3, 5), 1e-9, "tight"),
               ("Euler-7, tight constants", (2, 3, 5, 7), 1e-9, "tight"),
               ("Euler-11, tight constants", (2, 3, 5, 7, 11), 1e-9, "tight")]
    out = []
    for lam_t in (0.175, 0.170, 0.165):
        for name, P, slack, gmode in configs:
            best = None
            for y0 in (0.06, 0.08, 0.10, 0.12):
                t0 = lam_t - y0 * y0 / 2
                g = math.exp(0.02 * y0) * (1 + 1e-9) if gmode == "tight" else 1.03
                t1 = time.time()
                F = F_single(N0, t0, y0, P, slack, g)
                if best is None or F > best[0]:
                    best = (F, y0, t0, time.time() - t1)
            out.append(dict(Lambda=lam_t, config=name, F=best[0], y0=best[1], t0=best[2]))
            print(f"Lambda={lam_t}  {name:30s}  best F(N0) = {best[0]:+.5f}  (y0={best[1]}, {best[3]:.0f}s/eval)", flush=True)
    json.dump(out, open("results/push165.json", "w"), indent=1)

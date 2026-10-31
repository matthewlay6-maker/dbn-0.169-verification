"""Analytic monotonicity of the tail quantity B(N) (Polymath15 Lemma 'largen' applied twice), checked in arb.

B(N) = (max(T1, T2) + 1.006 N^{-y} max(T3, T4)) * log(N/M),  M = N0l (sums up to M are in A), u = log N:
  T1 = M^{1-s1} b_M,  T2 = N^{1-s1} b_N,  T3 = M^{1-s2} b_M,  T4 = N^{1-s2} b_N,
  s1 = (1+y)/2 + (t/2) u - c,  s2 = s1 - y - c',   b_N = exp(t u^2 / 4).
Each of  g_i(u) = log[ (prefactor) T_i * (u - log M) ]  has derivative
  (a) T1:            -(t/2) log M + 1/(u - log M)
  (b) T2:            (1-y)/2 + c - (t/2) u + 1/(u - log M)
  (c) N^{-y} T3:     -y - (t/2) log M + 1/(u - log M)
  (d) N^{-y} T4:     (1-y)/2 + c + c' - (t/2) u + 1/(u - log M)
All four are decreasing in u (for u > log M), so if each is < 0 at u = log N1, every term, hence each max,
hence B, is non-increasing on [N1, inf): B(N) <= B(N1) for all N >= N1.
"""
import json, sys
import flint
from flint import arb

flint.ctx.prec = 200
ROWS = {"0.18 (ours)": (0.17, 0.1414, 50_000_000), "0.1787854 (Gomila row)": (0.16125, 0.18727, 50_000_000),
        "0.176 (ours)": (0.171, 0.10, 100_000_000)}
M = 200_000
out = {}
for name, (t, y, N1) in ROWS.items():
    t, y = arb(repr(t)), arb(repr(y))
    c, c2 = arb("1e-20"), arb("1e-15")
    u = arb(N1).log(); lM = arb(M).log()
    inv = 1 / (u - lM)
    d = {"a": -(t / 2) * lM + inv,
         "b": (1 - y) / 2 + c - (t / 2) * u + inv,
         "c": -y - (t / 2) * lM + inv,
         "d": (1 - y) / 2 + c + c2 - (t / 2) * u + inv}
    ok = all(v < 0 for v in d.values())
    out[name] = dict(derivatives_upper={k: float(v.upper()) for k, v in d.items()}, all_negative=ok)
    print(f"{name}: derivative bounds at u = log N1: " + ", ".join(f"{k}<={float(v.upper()):.4f}" for k, v in d.items())
          + f"  -> B non-increasing on [N1, inf): {ok}")
json.dump(out, open("results/tail_monotone.json", "w"), indent=1)

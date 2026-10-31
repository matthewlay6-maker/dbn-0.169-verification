"""Assemble a certified cloud run: recombine block pieces in exact rational arithmetic and check barrier coverage.
For each block [Nm, Np]: F >= 1 - ca1_upper - sum(S_upper over pieces) - Z_upper, where every *_upper value is an
upward-rounded double from arb (so each is >= the true value) and the sum is done exactly with Fractions.
The pieces of a block must tile [2, 30 Np + 1) exactly.  Usage: python assemble_parts.py results/cert176 0.171"""
import glob, json, sys
from fractions import Fraction as Fr

d, t0 = sys.argv[1], float(sys.argv[2])
cfg = json.load(open("cloud/tasks.json"))
need = cfg["need"]
rows = [json.loads(l) for f in sorted(glob.glob(f"{d}/task-*.jsonl")) for l in open(f) if l.strip().startswith("{")]
parts = {}
for r in rows:
    if r.get("mode") == "block":
        parts.setdefault((r["Nminus"], r["Nplus"]), []).append(r)
res, bad = [], []
for (a, b) in [tuple(x) for x in cfg["blocks"]]:
    ps = sorted(parts.get((a, b), []), key=lambda r: r["slo"])
    Dp = ps[0].get("D", 30) if ps else 30
    tiles = len(ps) > 0 and all(p.get("D", 30) == Dp for p in ps) and ps[0]["slo"] == 2 and ps[-1]["shi"] == Dp * b + 1 and all(ps[i]["shi"] == ps[i + 1]["slo"] for i in range(len(ps) - 1))
    if not tiles:
        bad.append((a, b, "pieces do not tile the n-range"))
        continue
    F = 1 - Fr(ps[0]["ca1_upper"]) - sum(Fr(p["S_upper"]) for p in ps) - Fr(ps[0]["Z_upper"])
    res.append(dict(Nminus=a, Nplus=b, D=Dp, pieces=len(ps), F_lower=float(F), ok=F > Fr(need)))
bars = sorted([r for r in rows if r.get("mode") == "barrier" and "ta" in r], key=lambda r: r["ta"])
bfail = [r for r in rows if r.get("mode") == "barrier" and r.get("fail")]
cover = bool(bars) and bars[0]["ta"] == 0.0 and all(b["ok"] for b in bars) \
    and all(bars[i]["t_reached"] >= bars[i + 1]["ta"] for i in range(len(bars) - 1)) and bars[-1]["t_reached"] >= t0
items = [r for r in rows if "returncode" in r]
summ = dict(n_blocks=len(cfg["blocks"]), blocks_assembled=len(res), blocks_bad=bad,
            min_F_lower=min(r["F_lower"] for r in res) if res else None, need=need,
            all_blocks_ok=len(res) == len(cfg["blocks"]) and all(r["ok"] for r in res),
            weakest=sorted(res, key=lambda r: r["F_lower"])[:3],
            barrier_slices=len(bars), barrier_fail_rows=bfail, barrier_covers=cover,
            barrier_steps=sum(b["steps"] for b in bars), barrier_min=min(b["min_abs_f"] for b in bars) if bars else None,
            nonzero_returncodes=[i for i in items if i["returncode"] != 0], vcpu_seconds=sum(i["seconds"] for i in items))
json.dump(summ, open(f"{d}/summary.json", "w"), indent=1)
print(json.dumps({k: v for k, v in summ.items() if k not in ("weakest",)}, indent=1))
print("weakest blocks:", summ["weakest"])

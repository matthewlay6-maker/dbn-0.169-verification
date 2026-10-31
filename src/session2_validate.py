"""Session 2a: check Theorem 1.3 (Polymath15) against the certified integral H_t from ht.py.
Pass criterion per point: |H_t/B_t - f_t| <= e_A + e_B + e_{C,0} (rigorously, upper < bound)."""
import json, time, math
import flint
from flint import arb, acb
import ht, pm15

rows = []
t00 = time.time()
for x in (300, 1000):
    for t in ("0.1", "0.2", "0.4"):
        for y in ("0", "0.2", "0.9"):
            t0 = time.time()
            # H decays like e^{-pi x/8}: give the integral enough bits to resolve it
            flint.ctx.prec = int(math.pi * x / 8 * 1.4427) + 200
            H = ht.H(arb(t), acb(arb(x), arb(y)))
            P = pm15.parts(t, x, y)
            R = H / P["B"]
            f = pm15.f_t(P)
            eab, ec0 = pm15.err_bounds(P)
            diff = abs(R - f).upper()
            ok = bool(diff < eab + ec0)
            rows.append(dict(x=x, t=t, y=y, N=P["N"], ratio=R.str(10, radius=False), f=f.str(10, radius=False),
                             abs_diff=float(diff), bound=float(eab + ec0), eAB=float(eab), eC0=float(ec0),
                             ok=ok, prec=flint.ctx.prec, s=time.time() - t0))
            print(f"x={x} t={t} y={y} |H/B-f|={float(diff):.2e} bound={float(eab+ec0):.2e} {'ok' if ok else 'FAIL'} {time.time()-t0:.0f}s", flush=True)
json.dump(dict(rows=rows, runtime_s=time.time() - t00), open("results/session2_validate.json", "w"), indent=1)
print(f"Theorem 1.3 check: {sum(r['ok'] for r in rows)}/{len(rows)} pass, {time.time()-t00:.0f}s")

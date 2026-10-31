"""Session 1 sanity checks: (1) H_0 vs (1/8) xi; (2) sign changes of H_0 -> 2*gamma_k, k <= 20.
Also: real sign changes of H_t for t = 0.2 in the same window (exploratory, certified signs)."""
import json, time
import flint
from flint import arb, acb
import ht

out = {}
t0 = time.time()
flint.ctx.prec = 400
chk1 = []
for z in ("0", "10", "28.27", "50+3j", "100"):
    zz = acb(complex(z.replace("j", "j"))) if "j" in z else acb(arb(z))
    h = ht.H(0, zz)
    x = ht.xi_check(zz)
    d = h - x
    chk1.append(dict(z=z, H=h.real.str(30, radius=False), xi8=x.real.str(30, radius=False),
                     diff_bound=float(abs(d).upper()), overlap=bool(h.overlaps(x))))
out["check1"] = chk1

ref = {}
for line in open("results/zeros_200_hiprec.txt"):
    k, g = line.split()
    ref[int(k)] = arb(g)


def sign(t, z):
    v = ht.H(t, acb(z)).real
    if v > 0:
        return 1
    if v < 0:
        return -1
    return 0  # undetermined


def zeros_in(t, zmax, step):
    flint.ctx.prec = 160
    z = arb(0)
    s_prev = sign(t, z)
    found = []
    k = 1
    while float(z) < zmax:
        z2 = z + step
        s2 = sign(t, z2)
        if s_prev and s2 and s2 != s_prev:
            a, b, sa = z, z2, s_prev
            for _ in range(90):  # bisect to ~1e-27
                m = (a + b) / 2
                sm = sign(t, m)
                if sm == 0:
                    break
                if sm == sa:
                    a = m
                else:
                    b = m
            found.append((a, b))
        if s2:
            s_prev = s2
        z = z2
    return found


zs = zeros_in(0, 2 * 77.5, arb("0.2"))
chk2 = []
for k, (a, b) in enumerate(zs[:20], 1):
    g2 = 2 * ref[k]
    chk2.append(dict(k=k, lo=a.str(28, radius=False), hi=b.str(28, radius=False),
                     two_gamma=g2.str(28, radius=False),
                     contains=bool(float(a) <= float(g2) <= float(b)) or abs(float((a + b) / 2 - g2)) < 1e-20,
                     abs_err=abs(float((a + b) / 2 - g2))))
out["check2"] = dict(n_sign_changes=len(zs), rows=chk2)

zt = zeros_in(0.2, 2 * 77.5, arb("0.2"))
out["t0.2_real_sign_changes"] = dict(count=len(zt), first=[float((a + b) / 2) for a, b in zt[:20]])
out["runtime_s"] = time.time() - t0
json.dump(out, open("results/session1.json", "w"), indent=1)
ok1 = all(r["overlap"] and r["diff_bound"] < 1e-30 for r in chk1)
ok2 = len(chk2) == 20 and all(r["abs_err"] < 1e-20 for r in chk2)
print(f"check1 {'PASS' if ok1 else 'FAIL'} (max diff {max(r['diff_bound'] for r in chk1):.1e}); "
      f"check2 {'PASS' if ok2 else 'FAIL'} ({len(zs)} sign changes, max |mid-2g| {max(r['abs_err'] for r in chk2):.1e}); "
      f"t=0.2: {len(zt)} real sign changes; {time.time()-t0:.0f}s")

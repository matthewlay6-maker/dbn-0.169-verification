import math
import numpy as np


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



"""Polymath15 (arXiv:1904.12438) Theorem 1.3 in arb ball arithmetic.

H_t(x+iy)/B_t(x+iy) = f_t(x+iy) + O(e_A + e_B + e_{C,0})   on 0 < t <= 1/2, 0 <= y <= 1, x >= 200.
Equation numbers refer to the paper: (6) M_0, (9) alpha, (10) M_t, (11) B_t, (14)-(19) f_t,
gamma, s_*, kappa, N, and the explicit bounds (23) [e_A + e_B] and (24) [e_{C,0}].
"""
import flint
from flint import arb, acb

I = acb(0, 1)


def M0(s):
    pi = arb.pi()
    return (s * (s - 1) / 2) / 8 * (-(s / 2) * pi.log()).exp() * (2 * pi).sqrt() \
        * ((s / 2 - acb("0.5")) * (s / 2).log() - s / 2).exp()


def alpha(s):
    return 1 / (2 * s) + 1 / (s - 1) + (s / (2 * arb.pi())).log() / 2


def Mt(t, s):
    return (t / 4 * alpha(s) ** 2).exp() * M0(s)


def parts(t, x, y):
    t, x, y = arb(t), arb(x), arb(y)
    sp = (1 + y - I * x) / 2          # s_+ argument of B_t
    sm = (1 - y + I * x) / 2
    B = Mt(t, sp)
    gam = Mt(t, sm) / B
    sstar = sp + t / 2 * alpha(sp)
    kap = t / 2 * (alpha(sm) - alpha((1 + y + I * x) / 2))
    N = int(((x / (4 * arb.pi()) + t / 16).sqrt()).floor().unique_fmpz())
    return dict(t=t, x=x, y=y, B=B, gamma=gam, sstar=sstar, kappa=kap, N=N)


def f_t(P):
    t, y, N = P["t"], P["y"], P["N"]
    s1 = acb(0)
    s2 = acb(0)
    e1 = P["sstar"]
    e2 = P["sstar"].conjugate() + P["kappa"]
    for n in range(1, N + 1):
        ln = arb(n).log()
        b = (t / 4 * ln * ln).exp()
        s1 += b * (-e1 * ln).exp()
        s2 += b * (y * ln).exp() * (-e2 * ln).exp()
    return s1 + P["gamma"] * s2


def err_bounds(P):
    """Upper bounds (arb) for e_A + e_B, eq. (23), and e_{C,0}, eq. (24)."""
    t, x, y, N = P["t"], P["x"], P["y"], P["N"]
    pi = arb.pi()
    g = abs(P["gamma"])
    rs = P["sstar"].real
    kk = abs(P["kappa"])
    eab = arb(0)
    for n in range(1, N + 1):
        ln = arb(n).log()
        b = (t / 4 * ln * ln).exp()
        L2 = (x / (4 * pi * n * n)).log() ** 2
        eab += (1 + g * (kk * arb(N).log()).exp() * (y * ln).exp()) * b * (-rs * ln).exp() \
            * (((t * t / 16 * L2 + arb("0.626")) / (x - arb("6.66"))).exp() - 1)
    lx = (x / (4 * pi)).log()
    ec0 = (-(1 + y) / 4 * lx).exp() * (-t / 16 * lx * lx
                                         + arb("1.24") * (arb(3) ** y + arb(3) ** (-y)) / (N - arb("0.125"))
                                         + (3 * abs(acb(lx, pi / 2)) + arb("10.44")) / (x - 12)).exp()
    return eab.upper(), ec0.upper()

"""Certified evaluation of H_t(z) = int_0^inf e^{tu^2} Phi(u) cos(zu) du in arb.

Phi(u) = sum_{n>=1} (2 pi^2 n^4 e^{9u} - 3 pi n^2 e^{5u}) exp(-pi n^2 e^{4u}).
Every returned value is an acb ball that contains the true value:
  - Phi: truncated at n <= M plus a rigorous geometric bound on the tail n > M
    (valid where c = Re e^{4u} > 0; otherwise a non-finite ball, so the integrator subdivides);
  - integral: acb.integral (rigorous) on [0, U] plus an explicit bound on [U, inf).
"""
import math
import flint
from flint import arb, acb

PI = None


def _pi():
    return arb.pi()


def phi(u):
    """Phi(u) for an acb ball u, with rigorous truncation error."""
    pi = _pi()
    e4 = (4 * u).exp()
    c = e4.real.lower()
    if not (c > 0):
        return acb(arb("nan"), arb("nan"))
    prec = flint.ctx.prec
    cf = float(c)
    M = int(math.ceil(math.sqrt((prec + 30) * math.log(2) / (math.pi * cf)))) + 1
    if M > 400:
        return acb(arb("nan"), arb("nan"))
    e9 = (9 * u).exp()
    e5 = (5 * u).exp()
    s = acb(0)
    for n in range(1, M + 1):
        n2 = n * n
        s += (2 * pi ** 2 * n2 * n2 * e9 - 3 * pi * n2 * e5) * (-pi * n2 * e4).exp()
    # tail n > M: |term_n| <= (2 pi^2 n^4 |e9| + 3 pi n^2 |e5|) exp(-pi n^2 c),
    # successive ratio <= ((M+2)/(M+1))^4 exp(-pi (2M+3) c) =: q < 1
    n = M + 1
    first = (2 * pi ** 2 * n ** 4 * abs(e9).upper() + 3 * pi * n ** 2 * abs(e5).upper()) * (-pi * n * n * c).exp()
    q = (arb(M + 2) / (M + 1)) ** 4 * (-pi * (2 * M + 3) * c).exp()
    if not (q < 1):
        return acb(arb("nan"), arb("nan"))
    R = (first / (1 - q)).upper()
    return s + acb(arb(0, R), arb(0, R))


def _cutoff(prec):
    return math.log((prec + 40) * math.log(2) / math.pi) / 4 + 0.3


def tail_bound(t, z, U):
    """Rigorous bound on |int_U^inf e^{tu^2} Phi(u) cos(zu) du| for real t, U >= 1.
    For u >= U: |Phi(u)| <= 2 (2 pi^2 + 3 pi) e^{9u} e^{-pi e^{4u}} (n = 1 dominates; ratio
    of later terms < 1/2 since e^{4u} > 50), |cos(zu)| <= e^{|Im z| u}.  With
    E(u) = 9u + t u^2 + |Im z| u - pi e^{4u}, E'(u) <= -1 on [U, inf) once
    4 pi e^{4U} >= 10 + |Im z| + 2|t|U (E'' < 0 there), so the tail is <= C e^{E(U)}."""
    pi = _pi()
    U = arb(U)
    t = arb(t)
    y = abs(acb(z).imag).upper()
    assert float((4 * U).exp()) > 50
    assert 4 * pi * (4 * U).exp() >= 10 + y + 2 * abs(t) * U + 2 * abs(t)
    E = 9 * U + t * U * U + y * U - pi * (4 * U).exp()
    return (2 * (2 * pi ** 2 + 3 * pi) * E.exp()).upper()


def H(t, z, U=None):
    """H_t(z) as an acb ball (t real, z complex)."""
    prec = flint.ctx.prec
    if U is None:
        U = _cutoff(prec)
    t = arb(t)
    z = acb(z)
    f = lambda u, analytic: (t * u * u).exp() * phi(u) * (z * u).cos()
    val = acb.integral(f, 0, arb(U), rel_tol=arb(2) ** (-prec + 10))
    R = tail_bound(t, z, U)
    return val + acb(arb(0, R), arb(0, R))


def xi_check(z):
    """(1/8) xi(1/2 + iz/2) as an acb ball."""
    pi = _pi()
    s = acb(0.5) + acb(0, 1) * acb(z) / 2
    xi = s * (s - 1) / 2 * pi ** (-s / 2) * (s / 2).gamma() * s.zeta()
    return xi / 8

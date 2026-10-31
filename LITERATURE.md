# Literature: upper bounds on Λ (checked 2026-10-30)

| Bound | Source | Status | Notes |
|---|---|---|---|
| Λ ≥ 0 | Rodgers–Tao (2018/2020) | Published | Newman's conjecture; RH ⇔ Λ = 0 |
| Λ ≤ 0.22 | Polymath15, arXiv:1904.12438 (2019) | Published | Effective A+B−C approximation of H_t; barrier method; three x-regions |
| Λ ≤ 0.2 | Platt–Trudgian (2021) | **Published; best refereed unconditional bound** | Uses RH verified to T = 3,000,175,332,800 |
| Λ ≤ 0.1787854 | Preprint / blog writeup (2026) | Unrefereed claim | Computer-assisted |
| Λ ≤ 3/16 = 0.1875 | teorth/optimizationproblems (same author as below) | Self-archived certificate | Used as the anchor for the 0.16 / 0.1595 claims |
| Λ < 0.16, Λ < 0.1595 | teorth/optimizationproblems PRs #208, #210 | Merged into a records table 2026-10-18; **no formal review** | "Fully AI-derived" by the author's account. Collision exclusion on t ∈ [0.1595, 3/16]; interior x ≤ X and exterior N ≥ 690988 modules; X = 6,000,000,185,827; uses Platt–Trudgian at height X/2. A reviewer declined to accept it: the bundle shows the parts fit together, not that they are true; the anchor is the author's own result; cited "banks" are missing. |
| Λ ≤ 0.148612 | Planat, arXiv:2610.18214 (2026-10-13) | Preprint, **conditional** on stated "source, localisation and contact premises" | Rational Bernstein certificates, 192/256-bit enclosures, exact-rational auditor |
| Λ ≤ 0.13963 | arXiv:2610.12803 | Preprint, **conditional** on no zeta zeros with Re s > 7/8 (from an unrefereed preprint) | — |
| Λ_L ≤ 0.2978 | Freire et al., arXiv:2610.27343 | Preprint | Dirichlet L-function analogue; does not bound the classical Λ |

## Takeaways
- The refereed state of the art is still **Λ ≤ 0.2**.
- Several sub-0.2 claims appeared in Oct 2026. They are unrefereed, conditional, or
  certificate bundles nobody has independently re-derived.
- So an independent re-derivation of the key inequalities in our own code has real
  value right now, more than another attempt to lower the bound.

## Sources
- https://en.wikipedia.org/wiki/De_Bruijn%E2%80%93Newman_constant
- https://arxiv.org/abs/1904.12438 (Polymath15)
- https://terrytao.wordpress.com/2019/04/30/11075/
- https://github.com/teorth/optimizationproblems/pull/208
- https://github.com/teorth/optimizationproblems/pull/210
- https://arxiv.org/abs/2610.18214
- https://arxiv.org/html/2610.12803
- https://arxiv.org/abs/2610.27343
- https://www.judegomila.com/posts/riemann-lambda-0.1787854

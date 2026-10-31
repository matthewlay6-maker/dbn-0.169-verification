# Audit note: the claimed bound Λ < 319/2000 = 0.1595

Status: **independent read-only review plus numerical probes. Not a refutation and not a confirmation.**
Subject: the certificate merged as teorth/optimizationproblems PR #210, on top of #208 (0.16) and
the same author's earlier Λ ≤ 3/16. We read its `math/ASSEMBLY.md` (pinned commit `cc81078`) on the
web. Nothing from the bundle was downloaded or executed, and none of its tables were reused.

## 1. What the certificate claims (as we read it)
Collision exclusion. If Λ > 0, then H_Λ has a multiple real zero. Combined with the anchor Λ ≤ 3/16,
it suffices to show that H_t has no multiple real zero for t ∈ [0.1595, 0.1875]. The argument has
three parts:

| Part | Range | How it is supported in the bundle |
|---|---|---|
| Interior | \|x\| ≤ X = 6,000,000,185,827 | Platt–Trudgian RH to X/2, plus cited "separator/simplicity" results and boundary numerics (banks 590, 596, 599, 616, 617). Imported, not reproduced. |
| Exterior, finite | N = ⌊√(x/4π + t/16)⌋ ∈ [690988, 2299999], i.e. x up to ≈ 6.6×10¹³ | 42 numerical "leaves" with margins *asserted* to exceed 1/250000. Outputs not shown. |
| Exterior, infinite | N ≥ 2,299,999 | Analytic "supersolution" bound J(z) ≤ B∞. Monotonicity in N asserted. |
| Tail | \|x\| ≥ 4π e^{512000/319} | Cited results 588/590 |
| Anchor | Λ ≤ 3/16 | The same author's earlier certificate (not refereed) |

Key quantitative claim: at every real zero with x ≥ X and t ∈ [a, b], |u_x| > 1/400000, where
u = H_t/|M_t|.

## 2. Why the exterior module carries the whole result
Polymath15 §10 (Table 1, reproduced in our `cost_map.py`) and its heuristic argument say a bound
Λ ≤ c requires excluding RH violations up to a height that grows roughly like exp(C/c).
- Polymath15's own table puts **Λ ≤ 0.16 at X ≈ 2×10¹⁵**, about 330× beyond Platt–Trudgian's X = 6×10¹².
- A zeta zero off the line at height T with y ≈ 1 would move under the heat flow at speed ~log T. It would
  meet the real axis, as a collision, at time t ~ c/log T.
- So for T between ~10¹³ and ~10¹⁵, such a collision would fall inside t ∈ [0.1595, 0.1875] at x ≈ 2T,
  which is exactly the exterior region.

Hence the exterior claim, if true, rules out RH counterexamples of that kind at heights well beyond
any published verification. It does so with 42 finite leaves plus one analytic bound. That is not
impossible in principle: it is a statement about H_t at t ≈ 0.16, not about ζ. But it is exactly the
point where Polymath15 expected a "major breakthrough" to be needed. **It is the step that most needs
an independent check, and the bundle does not show the leaf margins or the analytic B∞ derivation.**

## 3. What we checked independently
| window start x | N | real zeros found | min \|f′\| at zeros | median | e_{C,0} bound |
|---|---|---|---|---|---|
| X = 6,000,000,185,827 | 690988 | 428 | 7.71 | 12.83 | 8.9e-07 |
| X + 10⁶ | 690988 | 428 | 8.03 | 12.95 | 8.9e-07 |
| X + 10⁸ | 690994 | 428 | 7.73 | 12.95 | 8.9e-07 |
| 2×10¹³ | 1261566 | 447 | 8.36 | 13.34 | 3.4e-07 |
| 4×10¹³ | 1784124 | 458 | 9.08 | 13.85 | 1.9e-07 |
| 6.6×10¹³ | 2291748 | 466 | 9.08 | 14.09 | 1.3e-07 |

t = 319/2000, windows of length 200, grid step 0.02, float64 with reduced phases. 2655 real zeros in total.
The smallest |d/dx (H_t/B_t)| at any zero is 7.71, about 3×10⁶ times the certificate's
claimed floor 1/400000. The sample is **consistent with the exterior claim**. It cannot confirm a statement for all x ≥ X.

Context from our cost map: the established Polymath15 method, even on optimistic assumptions, reaches 0.1595 only
with RH verified to ≈ 3×10¹³ (10× Platt–Trudgian). Our own full computation at Platt–Trudgian's height gives
0.18 (ball-arithmetic certification in progress), and 0.175 fails. So the claim requires an argument genuinely beyond that method, which is consistent with §2.

## 4. What would settle it
1. Publication of the 42 leaf outputs: the polydisc minima m(s), the D̄(s) values, and the leaf boundaries.
   Each leaf could then be re-derived independently. Our Theorem 1.3 / f_t implementation is ready for this.
2. A written derivation of the infinite-branch bound (D4) and its monotonicity in N, which the
   bundle asserts.
3. The cited banks (588, 590, 596, 599, 616, 617, 713, 729, 744, 745, 754, 755, 767, 771, 774, 779,
   786, 791) as checkable documents.
4. A re-check of the anchor Λ ≤ 3/16, which the chain needs.

## 5. Summary
The logic of collision exclusion is standard. The certificate's weight rests on an exterior module whose
key numbers are asserted, not shown. If that module holds, it implies a statement Polymath15 judged far
out of reach. Our probes found nothing inconsistent with the claim in the windows we sampled, and
sampling cannot confirm a "for all x" claim. Until items 1–4 are available, we would describe
Λ < 0.1595 as **unverified**, and Λ ≤ 0.2 (Platt–Trudgian) as the best independently checkable bound.

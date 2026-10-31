# de Bruijn–Newman: independent verification — results

All results are computed, or certified under the inputs listed. Nothing here is a new theorem about Λ,
and nothing bears on RH beyond "consistent with". Code: `ht.py`, `pm15.py`, `session1.py`,
`session2_validate.py`, `repro_pm15.py`, `barrier_pm15.py`, `cost_map.py`, `audit_spotcheck.py`.
Paper source: arXiv:1904.12438 (Polymath15), read from its LaTeX in `refs/src/`.

## Session 1 — certified H_t (arb ball arithmetic)
H_t(z) = ∫₀^∞ e^{tu²} Φ(u) cos(zu) du is evaluated with a rigorous truncation of the Φ series
(geometric tail bound) and a rigorous [U, ∞) bound. The integral uses `acb.integral` (rigorous).

**Check 1: H_0 vs (1/8)ξ(1/2 + iz/2)** (400 bits). Requirement: < 1e-30 with overlapping enclosures.

| z | H_0(z) | \|H_0 − ξ/8\| bound | enclosures overlap |
|---|---|---|---|
| 0 | 0.0621400972735392637 | 4.6e-120 | yes |
| 10 | 0.0344437496680255240 | 2.7e-120 | yes |
| 28.27 | −4.74989032118018463e-8 | 1.8e-120 | yes |
| 50+3i | −2.19377651835089208e-7 (real part) | 1.6e-120 | yes |
| 100 | 3.95274390744736142e-16 | 1.6e-120 | yes |

**PASS.**

**Check 2: real sign changes of H_0 on [0, 155]** (160 bits, certified signs, bisection). Found exactly
20 sign changes. Their midpoints match 2γ_k (k ≤ 20, mpmath at 28 digits) to ≤ 4.6e-23. **PASS.**
Exploratory: H_{0.2} also has 20 real sign changes on [0, 155]; the first is at 28.154 vs 2γ₁ = 28.270.
Runtime 435 s.

## Session 2 — Polymath15 Theorem 1.3 and an independent reproduction of Λ ≤ 0.22
### 2a. Effective approximation checked against the certified integral
At 18 points (x ∈ {300, 1000}, t ∈ {0.1, 0.2, 0.4}, y ∈ {0, 0.2, 0.9}), |H_t/B_t − f_t| ≤ e_A + e_B + e_{C,0}
held every time (18/18). Actual error / bound ranged from 0.24 to 0.53.

### 2b. Reproduction of Section 8 (t₀ = y₀ = 0.2, X = 6×10¹⁰ + 83952 − 0.5)

| Step | Paper | Ours | Status |
|---|---|---|---|
| (i) RH to X/2 ≈ 3×10¹⁰ | Platt | — | **external input** (now also covered by Platt–Trudgian 2021) |
| Sweep: e_A + e_B | ≤ 3.444e-9 | ≤ 1.36e-8 (arb) | holds; see anomaly B1 |
| Sweep: e_{C,0} | ≤ 1.249e-3 | ≤ 1.248358e-3 (arb) | reproduced |
| Sweep total | ≤ 1.25e-3 | ≤ 1.248371e-3 | **reproduced** |
| (c) A | ≤ 1.88 | **1.9063** (arb; float64 agrees) | **does not reproduce**; see anomaly B2 |
| (c) B at N₁ | ≤ 0.075 | 0.0740; decreasing on a log grid to 10⁴⁰ | reproduced (numerical monotonicity) |
| (c) conclusion | \|f_t − 1\| ≤ 0.955 | \|f_t − 1\| ≤ 0.9804, so \|f_t\| ≥ 0.0196 > 1.25e-3 | **holds after repair** |
| (b) right edge | uses 0.955 | with 0.9804: dist(E f, (−∞,0]) ≥ 0.0098 > 2.05e-3 | **holds after repair** |
| (b) F on [69098, 8e4] | 0.0263 | 0.02019 | holds (> 2.15e-3) |
| (b) F on [8e4, 1.1e5] | 0.0470 | 0.04082 | holds |
| (b) F on [1.1e5, 2.2e5] | 0.093 | 0.08673 | holds |
| (b) F on [2.2e5, 1.5e6] | 0.060 | 0.05255 | holds |
| (a) barrier, t ∈ [0, 0.2] | winding 0 | winding 0 in all 206 t-steps; min \|f_t\| = 1.55 on ∂R; mesh ≤ 14,746 points | **reproduced** |

Our claim (b) is stricter than the paper's. For every N in each block it uses the *true* truncated
mollified coefficients (n/d ≤ N), updated exactly as N grows, and takes the minimum over N. The
minimum is always at the right end of the block.

The barrier uses a 40-term Taylor expansion about the rectangle centre, with phases (X₀/2)·log n
reduced mod 2π in arb. Explicit allowances cover the Taylor tail, float rounding, and the α/κ variation.
It agrees with direct arb evaluation of f_t to ~1e-12.

**Conclusion of Session 2: Λ ≤ 0.22 reproduces independently, given Platt's RH verification to 3×10¹⁰,
after two minor repairs that leave every margin positive.** The (c) monotonicity of B and the
float64 parts of (b) and (a) carry explicit rounding allowances, but they are not ball arithmetic.

## Session 4 — computation toward Λ ≤ 0.18 (user-approved full run)
**Result: every check of Polymath15's Theorem 1.2 passes for t₀ = 0.17, y₀ = 0.1414, giving Λ ≤ t₀ + y₀²/2 =
0.17999698 < 0.18**, with RH verified to height 3,000,175,332,800 (Platt–Trudgian 2021) as the only external input.
Code: `full_run.py`, `barrier_gen.py`, `gonogo.py`, `gonogo2.py`. Results: `results/full_run.json`, `results/full_barrier.json`.

| Hypothesis / piece | Computed | Required |
|---|---|---|
| (i) no zeros with T ≤ X/2 | X = X₀ − 0.5 = 6000347387365.5, X/2 = 3,000,173,693,682.8 ≤ 3,000,175,332,800 | external (Platt–Trudgian) |
| Barrier location | X₀ = 6,000,347,387,366 (Euler product over p ≤ 11 ≈ 32.2), N = 691,008 constant on the barrier | — |
| Error budget, barrier (all t ∈ [0, t₀]) | e_A + e_B + e_{C,0} ≤ 4.650e-04 | — |
| Error budget, region (ii) (t = t₀), times \|E_{t,5}\| ≤ 1.6290 | 3.494e-07 | — |
| (iii) barrier [X₀ ± 0.5] × [y₀, 1], t ∈ [0, 0.17] | winding 0 in all 928 time steps; min \|f_t\| on ∂R = 1.4988; mesh ≤ 121,806 points | > 4.6e-04 + mesh/time slack |
| (ii) lower edge y = y₀, N ∈ [691,008, 50,000,000] | 110 blocks (ratio 1.04); min F = 0.1101 (first block) | ≥ 3.5e-07 |
| (ii) upper edge y = 1 | \|f − 1\| ≤ 0.610; dist(E f, (−∞,0]) ≥ 0.210 | > 3.5e-07 |
| (ii) left edge x = X₀ − 0.5 (t = t₀) | min \|f\| = 1.471, \|arg f\| ≤ 0.474; dist ≥ 0.789 | > 3.5e-07 |
| (ii) right edge N = N₁ and region N ≥ N₁ | A = 1.8252, B ≤ 0.0354: \|f − 1\| ≤ 0.8606, \|f\| ≥ 0.1394; right-edge dist ≥ 0.0749 | > 6.5e-11; > 1.1e-10 |

Go/no-go at the same location: Λ = 0.18 passes on two parameter routes (block F = 0.129 and 0.104). Λ = 0.175 **fails**
on both (F < 0), so the method's limit at this verification height is between 0.175 and 0.18.

**Status and caveats (read before quoting the number):**
- **Computed, not peer reviewed.** It depends on Polymath15's proven Theorem 1.2, Theorem 1.3 / Prop. "estimates",
  the derivative Lemma and the monotonicity Lemma (vi), as published, plus Platt–Trudgian's RH verification.
- The long sums (lower-edge blocks: up to 1.5×10⁹ terms; barrier moments) are float64 with explicit allowances:
  1e-6 absolute for blocks (worst-case relative rounding ~2e-7); Taylor tail, rounding and α/κ variation for the barrier.
  These are **not** ball arithmetic.
- Monotonicity of the tail quantity B in N was checked on a log grid from 5×10⁷ to ~10⁴⁰, not proved; the paper argues
  the analogous statement in a footnote.
- The memory-light block bound treats coefficients with N₋ < n/d ≤ N₊ in absolute value (worst case). It is valid for
  every N in the block, and stricter than the paper's untruncated coefficients.
- It is independent of, and weaker than, the unrefereed claims 0.1787854, 0.16 and 0.1595. Its value is transparency:
  every step is in this folder and re-runs in about 1.5 h on a desktop.


### Session 4b — ball-arithmetic certification of Λ ≤ 0.18 (Cloud Run, 2026-10-30)
`cloud/dbn.c` (FLINT 3 / arb) re-ran every long sum in ball arithmetic. Cloud Run Job `dbn-certify`, 10 × 1 vCPU,
13.8 vCPU-h, ≈ $1. Results: `results/cert018/`, summary `results/cert018/summary.json`.

| Piece | Certified (arb) | Float pipeline |
|---|---|---|
| Lower edge, 110 blocks, N ∈ [691008, 5×10⁷] | min F ≥ 0.1100596527 (all blocks present) | 0.1100586513 (allowance 1e-6) |
| Barrier, 4 slices chaining 0 → 0.1731527 ≥ 0.17 | 928 steps, min \|f\| ≥ 1.496509, winding 0 | 928 steps, 1.4988 |
| Upper edge y = 1 | S ≤ 0.6096023175, dist ≥ 0.2098398 | 0.6096023 / 0.2098398 |
| Tail region N ≥ N₁ = 5×10⁷, and right edge | \|f\| ≥ 0.1394 > 6.5e-11; dist ≥ 0.0749158 | same |
| Left edge (28,136 points, rigorous slack 7.1e-4) | min \|f\| = 1.4706036, dist ≥ 0.79007 | 1.4706036 / 0.789 |
| Error budgets e_barrier, need | 4.6495287e-4, 3.494031e-7 | identical |

**Tail monotonicity: proved** (`tail_monotone.py`, arb). Write u = log N. B(N) is a sum of maxima of four terms T_i(u)·(u − log M).
Each term's log-derivative is an explicit function decreasing in u. All four are ≤ −0.80 at u = log N₁ for every row we use
(0.18, 0.1787854, 0.176), so B is non-increasing on [N₁, ∞) and B ≤ B(N₁).
**Status: Λ ≤ 0.18 is certified in ball arithmetic, given the published inputs below.**
**Assumed:** Polymath15 Theorems 1.2 and 1.3 and its lemmas; Platt–Trudgian.


## Session 5 — certified Λ ≤ 0.176 (2026-10-30)
**t₀ = 0.171, y₀ = 0.10, Λ ≤ t₀ + y₀²/2 = 0.176**, barrier X₀ = 6,000,347,387,366. External input: RH to 3,000,175,332,800 (Platt–Trudgian).
Found by bisection (`bisect_wall.py`): our method's wall at this height lies between 0.17550 (fails at every y₀) and 0.17556.
Float pipeline (`full_run.py`, local) and **ball-arithmetic certification** (Cloud Run Job `dbn-certify-176`, 237 tasks at ~100
concurrent vCPUs, 25.2 vCPU-h, ≈ $2) agree.

| Piece | Certified (arb) | Float |
|---|---|---|
| Barrier, t ∈ [0, 0.171], 24 slices chained | 1188 steps, min \|f\| ≥ 1.439681, winding 0 | 1176 steps, 1.4253 |
| Lower edge, 127 blocks N ∈ [691008, 10⁸] (213 pieces, exact rational recombination) | min F ≥ 0.015049781 (first block) | 0.0150488 |
| Upper edge y = 1 | S ≤ 0.605148, dist ≥ 0.210277 | — |
| Left edge (29,493 points) | min \|f\| = 1.467262, dist ≥ 0.78099 | 1.4673 |
| Tail N ≥ 10⁸, right edge | \|f\| ≥ 0.0768 > 2.0e-11; dist ≥ 0.04036 | — |
| Error budgets | e_barrier = 6.141666e-04, need = 4.4383e-07 | same |
| Tail monotonicity of B | proved (`tail_monotone.py`) | — |

**Status: Λ ≤ 0.176 is certified in ball arithmetic, given Polymath15 Theorems 1.2 and 1.3, its lemmas, and Platt–Trudgian.
Not peer reviewed.**
Independent float verification of Gomila's row (Λ ≤ 0.17878502645 with y₀ = 0.18727 ≤ √(87677/2.5e6)): all six checks pass
(`results/g_run.json`); not yet ball-certified.


## Session 6 — certified Λ ≤ 0.173 (Euler-11 near N₀, split rectangle; 2026-10-31)
**t₀ = 0.168, y₀ = 0.10, Λ ≤ 0.173**, barrier X₀ = 6,000,347,387,366, N₁ = 2×10⁸. External input: Platt–Trudgian.
New rigorous ingredients:
1. Mollifier E_{t,P} with P = {2,3,5,7,11} on region A (N ∈ [691008, 932860]) and P = {2,3,5} on region B (N ∈ [932860, 2×10⁸]).
2. Tight constants from Polymath15 Thm 1.3 (20)–(21): σ slack 10⁻⁹ instead of 0.001, and \|γ\|nʸ ≤ e^{0.02y₀}(1+10⁻⁹)(n/N)ʸ instead of 1.03.
3. The argument principle applied separately on the two sub-rectangles. The internal vertical edge at x_s = 10,935,616,520,035
   is certified with both mollifiers.
Cloud: `dbn-certify-173`, 480 tasks at ~100 vCPUs, 45.9 vCPU-h (≈ $3).

| Piece | Certified (arb) |
|---|---|
| Barrier, t ∈ [0, 0.168], 24 slices | 1188 steps, min \|f\| ≥ 1.439681, winding 0 |
| Lower edge, region A (Euler-11, 31 blocks) | min F ≥ 0.024800 (first block); need 5.47e-7 |
| Lower edge, region B (Euler-5 tight, 137 blocks) | all positive (assembled, all 168 blocks tile exactly) |
| Upper edge y = 1 | S ≤ 0.6259935; dist ≥ 0.18358 (A), 0.19662 (B) |
| Left edge (A) | min \|f\| = 1.47737, dist ≥ 0.72478 |
| Internal vertical edge x_s | min \|f\| = 0.97779; dist ≥ 0.48941 (A), 0.52230 (B) |
| Tail N ≥ 2×10⁸ and right edge (B) | \|f\| ≥ 0.1413 > 6.1e-12; dist ≥ 0.07428; B monotone (proved) |

**Status: Λ ≤ 0.173 is certified in ball arithmetic, given the published inputs. Not peer reviewed.**
Scouting results: multi-tap Euler dampers add only ≈ +0.011 (per-prime gains follow p^{−σ}; `results/stage0`);
a b-weighted long mollifier (y = 100) gives F > 0 already at Λ = 0.170 (cloud grid running: `levinson/`).

## Session 3B — cost map (Polymath15 method)
| X (RH verified to X/2) | optimistic Λ (criterion A) | with 0.02 margin (B) | Polymath15 Table 1 (verified) |
|---|---|---|---|
| 6e+10 | 0.2148 (t₀=0.214, y₀=0.04) | 0.2160 | — |
| 2e+11 | 0.2026 (t₀=0.202, y₀=0.04) | 0.2038 | — |
| 6e+11 | 0.1926 (t₀=0.188, y₀=0.1) | 0.1939 | — |
| 2e+12 | 0.1826 (t₀=0.182, y₀=0.04) | 0.1836 | 0.21 |
| 6e+12 | 0.1742 (t₀=0.173, y₀=0.04) | 0.1752 | — |
| 2e+13 | 0.1662 (t₀=0.161, y₀=0.1) | 0.1669 | 0.19 |
| 6e+13 | 0.1591 (t₀=0.157, y₀=0.06) | 0.1598 | 0.18 |
| 2e+14 | 0.1520 (t₀=0.150, y₀=0.06) | 0.1528 | — |
| 6e+14 | 0.1462 (t₀=0.144, y₀=0.06) | 0.1468 | — |

Criterion A checks only the hardest single point N = N₀ (uniform in T, Euler-5 mollified) against the bare error budget.
It is an **optimistic proxy**: it skips the full block coverage, the barrier, the tail and the safety margins.
It runs 0.02–0.03 below Polymath15's verified table at the same X (e.g. 0.183 vs 0.21 at 2×10¹²).
Platt–Trudgian height (X ≈ 6×10¹²): proxy 0.174. Our full computation (Session 4) gives 0.18 there (float64 with allowances; ball-arithmetic certification in progress), and 0.175 fails.
Proxy Λ at the claimed 0.1595 needs X ≈ 6×10¹³, i.e. ~10× the verified height.

## Session 3A — audit of Λ < 0.1595
See [AUDIT.md](AUDIT.md). Spot-check numbers:
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

## Anomalies
- **B1.** In the Polymath15 sweep proposition, the chain to e_A + e_B ≤ 3.444e-9 contains slips: a
  factor 1.023 where the paper's own steps give (1 + 1.021)·1.001 ≈ 2.023, and exponents 0.5999/0.80005
  where 0.4001/0.79995 follow. Recomputed: ≤ 1.36e-8. Still negligible against e_{C,0}.
- **B2.** Claim (c): A = 1.9063, not ≤ 1.88, so the stated |f_t − 1| ≤ 0.955 does not follow. The weaker
  0.9804 still gives both uses of (c) with margin. **The theorem Λ ≤ 0.22 is unaffected.**
- **B3.** Claim (b) F values: ours are ~0.006–0.008 below the paper's stated values on all four blocks.
  We did not find the cause; all are still ≥ 9× the threshold. The paper's β_n ignore the truncation
  n/d ≤ N. In our exact treatment that turned out harmless here, but the text does not justify it.
- **B4.** Barrier algorithm: the paper's time-step condition writes D_t|t′ − t|/(2n); we used the
  conservative D_t|t′ − t|.
- **B5.** The cost-map optimum sits at the smallest y₀ in our grid (0.04). The frontier is flat in y₀
  (Λ changes < 0.001 between y₀ = 0.04 and 0.10 at X = 6×10¹⁰).

## Stage 1, Step 1 — reproduction of the 0.1787854 pipeline (2026-10-31)
Repository `judegomila/dbn-lambda-01787854-candidate-audit`, commit **a74738deb6d5e0f76887cb36901da08b68dca705**,
cloned read-only to `ext/gomila-audit/`. Notes with file and line references: `stage1/pipeline_notes.md`.
Their `src/lemma_sweep_p235711.c` was compiled with their exact flags (`gcc -O3 -DTRIANGLE_WEIGHT`, FLINT 3.0.1, WSL) and run
with their exact shard boundaries and arguments (`6 220 14 0.005 t`, t = 16125/100000, y₀² = 350708/10⁷) on
N = 690988 … 692987.
**Result: 2,000 / 2,000 per-window floors match their logged `L12` values exactly (string-identical), 0 mismatches.**
Minimum 7.91366×10⁻⁷ at N = 690988, matching `CANDIDATE_PARAMETERS.md`. Runtime about 855 s per shard here vs about 80 s in
their logs (≈ 11× slower: older FLINT and WSL). Comparison: `stage1/repro_compare.json`.
Barrier reproduction: pending.

## Session 7 — certified Λ ≤ 0.169 (long mollifier; 2026-10-31)
**t₀ = 0.1658, y₀ = 0.08, Λ ≤ t₀ + y₀²/2 = 0.169**, barrier X₀ = 6,000,347,387,366, N₁ = 5×10⁸. External input: Platt–Trudgian.
New ingredients relative to Session 6:
- **Region A (N ∈ [691008, 1243814]):** b-weighted Selberg/Levinson-type long mollifier, M(s) = Σ_{d≤299, squarefree} λ_d d^{−s},
  λ_d = μ(d)·Π_{p|d} b_p·P(log(300/d)/log 300), P(x) = x − 1.805(x² − x). Frozen as exact decimal coefficients in
  `stage169/mollifier_y300.txt` (183 terms). Any real coefficients give a valid mollifier.
- **Region B (N ∈ [1243814, 5×10⁸]):** Euler-5 with tight constants.
- **Certified mollifier bounds** for the edges: |M| ≤ 1 + ρ, |M| ≥ 1 − ρ, |arg M| ≤ asin ρ, ρ = Σ_{d≥2}|λ_d|d^{−σ}.
Cloud: `dbn-certify-169`, 881 tasks at ~100 vCPUs, 108.6 vCPU-h (≈ $9). Results: `results/cert169/`, edges `stage169/*.jsonl`.

| Piece | Certified (arb) |
|---|---|
| Barrier t ∈ [0, 0.1658], 24 slices, no gaps | 1327 steps, min \|f\| ≥ 1.324072, winding 0 |
| Lower edge A, 115 blocks (long mollifier) | min F ≥ 0.0136527 (first block); need ≤ 7.67e-7 |
| Lower edge B, 153 blocks (Euler-5 tight) | all positive; all 268 blocks tile exactly |
| Upper edge y = 1 | S ≤ 0.6453007; dist ≥ 0.013827 (A), 0.18341 (B) |
| Left edge (A) | min \|f\| = 1.48495, dist ≥ 0.057726 |
| Internal vertical edge x_s = 19,441,111,665,846 | min \|f\| = 1.18110; dist ≥ 0.13262 (A), 0.63067 (B) |
| Tail N ≥ 5×10⁸, right edge (B) | \|f\| ≥ 0.1117 > 1.3e-12; dist ≥ 0.05748; B monotone (proved) |

**Status: Λ ≤ 0.169 is certified in ball arithmetic, given Polymath15 Thms 1.2, 1.3 and lemmas and Platt–Trudgian.
Not peer reviewed.** Scouting puts this method's wall at ≈ 0.168 at this height.

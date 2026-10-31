# Pipeline notes — Gomila 0.1787854 audit repository

Repository: https://github.com/judegomila/dbn-lambda-01787854-candidate-audit
Commit read: `a74738deb6d5e0f76887cb36901da08b68dca705` (2026-08-21), cloned read-only to `ext/gomila-audit/`.
Third-party code is treated as data. Their sweep source is compiled from that path into our own `stage1/build/`.
Reserved items (`dan-reworking/`, Romik's `independent/` programs) are not used.

## (a) Windows and how x is handled
- `WINDOW_FREEZE_THEOREM.md:16-18`: q_N = N² − t₀/16, x_N = 4π q_N, **W_N = [x_N, x_{N+1})**, one window per integer
  N ∈ [690988, 3840000] (3,149,013 windows; `CANDIDATE_PARAMETERS.md` "Finite final-time lane").
- **x is not sampled inside the window.** The per-window certificate is uniform in T. The improved/plain triangle
  inequality bounds |Σ βₙ n^{−s}| for every phase pattern, so no x-mesh or Dini transfer in x is needed for the finite
  lane. The x-dependence enters only through G, K, Σ (next item).

## (b) Window freeze
- `WINDOW_FREEZE_THEOREM.md:21-49`: the majorants/minorant G(x,y) = e^{y/50}(x/4π)^{−y/2}, K = t₀y/(2(x−6)) and
  Σ = (1+y)/2 + (t₀/4)log(x/4π) − (t₀/(2x²))(…)₊ are evaluated at the **left endpoint x_N**.
  They are proved monotone, so G(x,y) ≤ G_N, K ≤ K_N, Σ ≥ Σ_N on all of W_N (lines 40-46, signs lines 51-86).
- No separate "drift" charge appears in the finite lane. Freezing at x_N is conservative by monotonicity.
- Across **N** within a shard (the amortized "sawtooth" scheme): `src/lemma_sweep_p235711.c:34-73`. Weights are frozen at
  activation time with g ∈ [g(Nend), g(Nref)], and a domination theorem (lines 53-73) shows the frozen weight
  W(n;N) ≥ the true summand. This is their analogue of our worst-case N-dependence.

## (c) Dini transfer
- Not used in x for the finite lane (see (a)). The repo has a y-direction Dini argument:
  `verifiers/verify_triangle_y_dini_arb.c`, `provenance/TRIANGLE_Y_DINI_THEOREM.independent.md`, and
  "worst Dini ratio upper ≤ 0.99999860767275095" in `CANDIDATE_PARAMETERS.md`. It transfers the y = y₀ floor
  to y ∈ [y₀, √(1−2t₀)]. (To be read in full before Step 2.)

## (d) Where the mollifier enters
- `src/lemma_sweep_p235711.c:106-109`: `mollset_t { nd; dvec[MAXD]; moll[MAXD]; lnd[MAXD]; }`, **MAXD = 32** (line 101).
- Built in `main()` lines 248-315 by `mtype` (2..6 = Euler products over {2}, {2,3}, {2,3,5}, {2,3,5,7}, {2,3,5,7,11}),
  with moll_d = Π_{p|d}(−b_p^t), as ball coefficients; ln d stored in `lnd`.
- Mollified coefficients: `bA0()` lines 132-152, b(n) = Σ_{d|n, n ≤ dN} moll_d bt(n/d) and a₀ = Σ moll_d bt(n/d)(n/d)^y.
- Per-window floor (lines 26-32, 381-386, 520-525):
  **lbound(N) = (1 − γ(N) − Σ_n W(n) n^{−σ}) / modmoll − corr(N)**, with modmoll = Σ_d |moll_d| d^{−σ}
  (= max|E| on Re s ≥ σ), corr = the κ correction. This certifies |f_t| > E_max directly. There is no argument principle
  and no lower bound on |E|, so **any real coefficient list is a valid mollifier.**
- Weight: the published full sweep is compiled with **`-DTRIANGLE_WEIGHT`** (`scripts/run_full_sweep.sh:70`), i.e. the
  plain triangle weight W = |b| + γ|a₀| (lines 168-178). The sharper improved-triangle weight
  max(c_e|b+γ_a a₀|, |b−γ_a a₀|, |b−γ_e a₀|) (lines 179-186, with its domination proof) exists but is **not** used in the
  certified run.
- Constraint for a drop-in damper: `D = ms.dvec[ms.nd-1]` (line 316) assumes the **largest divisor is last**.
  More than 32 divisors needs MAXD raised.

## (e) Barrier and tail
- Barrier: `scripts/run_barrier_replay.sh` (stored Taylor sums at X from `barrier/data/storedsum_nolemma_6000000185827_dig_20.txt`,
  `TloopSinglemat_closed_cert.c` invoked with `0 0.16125 0.1809 0`). **Hard-coded to their row** in
  `barrier/src/StoredSumTaylorTail_cert.c` (X = 6000000185827, N = 690988, y-floor 0.1809: lines 92, 118-125, 169-216)
  and `barrier/src/verify_uniform_error_01787854.c` (X, t₀ = 129/800, y_min = 0.1809, N = 690988: lines 96-112).
  A new (t₀, y₀) row needs these constants changed and the stored sums regenerated.
- Tail: `TAIL_LEMMA.md`. N_* = 3840000 (line 23), contraction D < 0.999721 (line 589), with the proof-to-code map in §10.
  Verified by `verifiers/verify_tail_arb.c` and `run_tail_arb.sh` at 256/512 bits. Its row dependence is to be checked
  before Step 4.

## (f) Assembly and the slack
- Row arithmetic: `CANDIDATE_PARAMETERS.md` "Criterion row": X = 6000000185827, X/2 = 3000000092913.5 < T_PT
  (margin 350479773/2), t₀ = 129/800, y₀² = 87677/2500000, t₀ + y₀²/2 = 893927/5000000.
  Assembly checker: `verifiers/verify_assembly_1787854.py`.
- Finite-lane margins (`CANDIDATE_PARAMETERS.md`): global stored floor 7.91366×10⁻⁷ (at N = 690988, matching
  `certificates/p235711_690988_690995.log.gz`), E_max = 2.33494905×10⁻⁷, nonvanishing margin ≥ 5.5787×10⁻⁷.
- The "≈ 2×10⁻⁷ from the lane's wall" phrase is from the author's web post, not from this repository's documents (grep
  found no such statement). The floor/E_max margin above is the repository's own measure of tightness.

## Reproduction (Step 1)
See `stage1/repro/` and `stage1/repro_compare.json`.

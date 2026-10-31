# dbn-0.169-verification

A certified computation of **Λ ≤ 0.169** for the de Bruijn–Newman constant, plus independent checks of earlier bounds.

> Subject to the correctness of Polymath15's Theorems 1.2 and 1.3 and lemmas
> ([arXiv:1904.12438](https://arxiv.org/abs/1904.12438)) and the Platt–Trudgian verification of RH to height
> 3,000,175,332,800, the Polymath criterion yields **Λ ≤ 0.169** (t₀ = 0.1658, y₀ = 0.08).
> Every certified inequality is decided in arb ball arithmetic. **Not peer reviewed.**

Write-up: [`paper/note.pdf`](paper/note.pdf) (source `paper/note.tex`). Full record of every run: [`results.md`](results.md).

## What is here
| Result | Status |
|---|---|
| Λ ≤ 0.169 | certified (ball arithmetic): barrier, all 268 lower-edge blocks, all edges, tail |
| Λ ≤ 0.173, 0.176, 0.18 | earlier certified rows (superseded) |
| Polymath15's Λ ≤ 0.22 | reproduced independently; errata in `paper/note.pdf` §4 |
| Gomila's Λ ≤ 0.1787854 | first 2,000 finite windows reproduced digit for digit (`stage1/`) |
| Method's limit at this height | ≈ 0.168 (float scouting) |

## The method in one paragraph
Polymath15's criterion needs a lower bound, uniform in T, for a mollified Dirichlet polynomial along the canopy's lower
edge. Earlier work used Euler-product mollifiers over a few small primes. Here the hardest range
(N ∈ [691008, 1243814]) uses a **b-weighted Selberg/Levinson-type long mollifier** over all squarefree d ≤ 299
(`stage169/mollifier_y300.txt`, 183 exact coefficients). It is combined with Polymath15's **improved triangle
inequality** and with constants taken directly from Theorem 1.3. An Euler-5 mollifier covers the rest, and a certified
internal edge joins the two regions.

## Layout
| Path | Contents |
|---|---|
| `cloud/dbn.c` | **Certified kernels (FLINT/arb):** block bounds with any mollifier (Euler sets, or a coefficient file via `DBN_MOLL`), barrier stepping, upper, right and tail edges, left and internal edges |
| `cloud/tasks.json`, `cloud/run.py`, `cloud/Dockerfile` | Exact work split and Cloud Run entrypoint for the 0.169 run (881 tasks) |
| `results/cert169/` | All task logs of the certified 0.169 run; `summary.json` |
| `stage169/` | Mollifier coefficients, edge certificates (`edges_*.jsonl`, `left_A`, `vert_*`), error budgets and tail |
| `results/cert173/`, `cert176/`, `cert018/` | Earlier certified rows |
| `src/` | Rigorous H_t (`ht.py`), Theorem 1.3 (`pm15.py`), Polymath15 reproduction, tail monotonicity proof, assembly, float scouting |
| `stage1/` | Reproduction of the 0.1787854 pipeline (notes with file and line references, comparison) |
| `AUDIT.md`, `LITERATURE.md` | Notes on the Λ < 0.1595 certificate; literature status as of 2026-10-31 |

## Reproduce
```
gcc -O2 -o dbn cloud/dbn.c -lflint -lmpfr -lgmp -lm                 # FLINT >= 3.0
# hardest block, long mollifier, tight constants (≈ 4 min):
DBN_MOLL=stage169/mollifier_y300.txt DBN_TIGHT=1 ./dbn block 691008 692390 0.1658 0.08
# edges for region A (long mollifier) and region B (Euler-5):
DBN_MOLL=stage169/mollifier_y300.txt ./dbn edges 6000347387366 0.08 0.1658 500000000 0.88825284
DBN_PRIMES=2,3,5 ./dbn edges 6000347387366 0.08 0.1658 500000000 0.88825284
# recombine the cloud pieces exactly:
PYTHONPATH=src python src/assemble_parts.py results/cert169 0.1658
```
The full certified run is about 109 vCPU-hours. `cloud/tasks.json` lists every task.

## Acknowledgement
The code and computations were developed with the assistance of an AI system (Claude, Anthropic). All results were
checked as described in `results.md`.

MIT License (code). © 2026 Matthew Davis.

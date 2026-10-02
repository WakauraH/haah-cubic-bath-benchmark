# Channel-resolved memory time of Haah's cubic code — reproduction package

Code, raw results and figure scripts for

> H. Wakaura, *Channel-resolved memory time of Haah's cubic code under reversible rate filtering and temporally correlated baths* (2026).

Everything quoted in the paper traces to a file in `results/`; `results/MANIFEST.md` records how each file was produced (script, seeds, date, and any later regeneration).

## Quick start

```bash
pip install -r requirements.txt      # numpy, scipy, numba, matplotlib, pytest
make check                           # ~15 s on 16 cores: construction, decoder, derived numbers, figures, two bit-reproducible MC points
make check-all                       # + L=2 exact master equation vs Gillespie, OU sigma->0 self-test (~8 min)
make numbers                         # regenerate results/derived_numbers.json
make figs                            # regenerate paper/figs/*.pdf
```

Tested with Python 3.13 on macOS (numpy 2.2/2.5, scipy 1.17/1.18, numba 0.61/0.68, matplotlib 3.10/3.11). The Monte Carlo kernels use `multiprocessing` with the `fork` start method, so Linux and macOS are supported; Windows is not.

## What `make check` verifies

| Test | What it recomputes | Compared against |
|---|---|---|
| code construction | Code 1 from the Bravyi–Haah corner table; all stabilizers commute; `k(L)` for `L=2..7` | `k(L)` formula of Haah 2011 Table II |
| decoder self-test | RG decoder on `L=5` (all weight-1 errors) and `L=9` (200 random weight-2/3 errors) | success thresholds in `scripts/haah_tmem.py` |
| derived numbers | `scripts/protocol_f.py` from the shipped `results/*.json` | `results/derived_numbers.json`, every number to 1e-9 |
| figures | `paper/figs/make_figs.py` | the five PDFs exist and are rebuilt |
| MC, flat bath | `L=5`, `β=4.5`, `n=100`, seeds `1000L+i` | `results/tmem_beta4.5_fast_flat.json` key `"5"`, mean and SEM to 1e-9 |
| MC, filtered bath | same with `(s4,s2,s0)=(1,1,0.001)` | `results/tmem_beta4.5_fast_s0_0.001_allL.json` key `"5"` |
| two-defect pseudomode | `scripts/pseudomode_twodefect.py` | `results/gap/B_twodefect_meanforce.json` to 1e-6 |

Bit-for-bit agreement of the Monte Carlo points holds because every run is seeded per sample (`seed = seed_offset + 1000*L + i`) and the rejection-free kernel is deterministic given the seed. Larger sizes reproduce the same way; see the timings below.

## Reproducing every result file

All commands run from `scripts/`. The output file name is `results/tmem_beta<beta><tag>.json` and **is overwritten for a given tag**, so use a fresh tag.

| Result | Command | Wall time (16 cores) |
|---|---|---|
| Code check | `python3 haah_check.py` / `python3 haah_code.py` | seconds |
| Decoder self-test | `python3 haah_tmem.py selftest` | seconds |
| L=2 exact master equation vs Gillespie | `python3 haah_davies.py` → `results/tier1_validation.json` | minutes |
| Flat bath, β=4.5, L=5–13 | `python3 haah_tmem_fast.py 4.5 5,7,9,11,13 100 _fast_flat` | 7 min |
| Flat bath beyond L\* | `python3 haah_tmem_fast.py 4.5 17 100 _fast_flat_L17`; `... 19 100 _fast_flat_L19` | 17 min; 24 min |
| Flat bath, β=5.0, L=5–11 | `python3 haah_tmem_fast.py 5.0 5,7,9,11 60 _fast_flat` | 5 min |
| Channel filter `(s4,s2,s0)` | `python3 haah_tmem_fast.py 4.5 5,7,9,11,13 100 _fast_s0_0.01 1,1,0.01` (and likewise for each row of Table I) | 2–8 min per series; 38 min for `s0=10` |
| Filtered bath beyond L\* | `python3 haah_tmem_fast.py 4.5 17 100 _fast_s0_0.01_L17 1,1,0.01`; `... 19 ... _L19 ...` | 14 min; 33 min |
| Creation-limited crossover | `python3 haah_tmem_fast.py 4.5 9 100 _gapC_s4_3e-4 3e-4,1,1` (and `1e-4`) | minutes |
| Effective temperature β′ | `python3 haah_tmem_fast.py 4.7 5,7,9,11,13 100 _fast_flat_lamb` (also 4.875, 5.06) | 13, 28, 57 min |
| OU-modulated diffusion rate | `python3 haah_tmem_ou.py 4.5 9 100 2.0 1,100,1e6` → `results/h2ou_beta4.5_L9_sigma2.0.json` (merges per τ_c) | 1–5 min per τ_c; the `τ_c=0.01` point of the σ=1 series took 5 h |
| Quenched OU points | `python3 haah_tmem_ou.py 4.5 9 200 2.0 1e10` (also L=7, 11 at σ=2; L=9 at σ=1) | 1–7 min |
| OU σ→0 self-test | `python3 haah_tmem_ou.py selftest` | minutes |
| Pseudomode models | `python3 pseudomode_gap.py`, `pseudomode_twodefect.py`, `pseudomode_cutoff_sweep.py`, `pseudomode_uniqueness.py` | seconds each |
| Derived numbers (Tables I–III, Sec. III) | `python3 protocol_f.py` → `results/derived_numbers.json` | seconds |
| Figures | `python3 ../paper/figs/make_figs.py` (from the repository root: `python3 paper/figs/make_figs.py`) | seconds |

The checks in `tests/` write `results/*_pytest_*.json`, which are ignored by git.

CLI of `haah_tmem_fast.py`: `<beta> <L,..> <nsamp> [tag] [s4,s2,s0] [seed_offset] [Ltqo|auto]`. `Ltqo` defaults to `L-2`; a trial decoding runs every `T_ec = e^{4β}/100`; the memory time of a run is the first trial at which the RG decoder fails or returns a logical error.

CLI of `haah_tmem_ou.py`: `<beta> <L> <nsamp> <sigma> <tau_c,..>`; seeds `7000+1000L+i`; `ξ(0)` is drawn from the seeded generator. Which OU records predate the seeding of `ξ(0)` (statistically valid, not bit-reproducible) is listed in `results/MANIFEST.md`.

## Layout

```
scripts/      haah_code.py (Code 1), haah_tmem.py (RG decoder, pure-Python kinetics), haah_tmem_fast.py (numba BKL kernel),
              haah_tmem_ou.py (OU-modulated rate), haah_davies.py (L=2 exact solution), pseudomode_*.py, protocol_f.py (all derived numbers)
results/      every result file quoted in the paper, run logs, MANIFEST.md (provenance), derived_numbers.json
paper/figs/   make_figs.py and the five figure PDFs
tests/        the checks run by `make check`
```

Files named `*fastcheck*`, `*n20*`, `*n40*`, `*LtqoLm2*`, `tmem_beta4.0.json` and `tmem_beta4.5_reviewA_flat.json` are validation runs quoted only in Appendix A; `MANIFEST.md` says which is which.

## License

MIT (see `LICENSE`).

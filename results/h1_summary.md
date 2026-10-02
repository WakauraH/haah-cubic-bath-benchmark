# H1 summary (β=4.5, n=100, L_tqo=L−2, Metropolis × channel filter)

| spec | L=5 | L=7 | L=9 | L=11 | L=13 | T(13)/T(5) |
|---|---|---|---|---|---|---|
| flat (1.0, 1.0, 1.0) | 2.28e+07 | 7.09e+07 | 1.20e+08 | 2.02e+08 | 2.26e+08 | 9.9 |
| s4_0.01 (0.01, 1.0, 1.0) | 4.62e+07 | 8.24e+07 | 1.14e+08 | 2.08e+08 | 2.66e+08 | 5.8 |
| s4_0.1 (0.1, 1.0, 1.0) | 2.22e+07 | 6.81e+07 | 9.61e+07 | 1.92e+08 | 2.46e+08 | 11.1 |
| s4_10 (10.0, 1.0, 1.0) | 2.22e+07 | 7.42e+07 | 1.16e+08 | 2.19e+08 | 2.77e+08 | 12.5 |
| s2_0.1 (1.0, 0.1, 1.0) | 6.85e+07 | 1.36e+08 | 2.62e+08 | 4.95e+08 | 5.49e+08 | 8.0 |
| s2_10 (1.0, 10.0, 1.0) | 1.03e+07 | 4.37e+07 | 6.85e+07 | 1.16e+08 | 1.63e+08 | 15.8 |
| s0_0.1 (1.0, 1.0, 0.1) | 9.78e+07 | 3.83e+08 | 5.37e+08 | 8.74e+08 | 1.26e+09 | 12.9 |
| s0_10 (1.0, 1.0, 10.0) | 9.34e+06 | 1.87e+07 | 3.45e+07 | 6.45e+07 | 6.06e+07 | 6.5 |

## Channel sensitivity  d ln T_mem / d ln s  (from s=0.1 and s=10 vs flat)

- s4 (|ΔE|=4 pair creation), L=9: slope = +0.01  (points s=[0.01, 0.1, 1.0, 10.0])
- s4 (|ΔE|=4 pair creation), L=13: slope = +0.00  (points s=[0.01, 0.1, 1.0, 10.0])
- s2 (|ΔE|=2 barrier/hook), L=9: slope = -0.29  (points s=[0.1, 1.0, 10.0])
- s2 (|ΔE|=2 barrier/hook), L=13: slope = -0.26  (points s=[0.1, 1.0, 10.0])
- s0 (ΔE=0 diffusion), L=9: slope = -0.60  (points s=[0.1, 1.0, 10.0])
- s0 (ΔE=0 diffusion), L=13: slope = -0.66  (points s=[0.1, 1.0, 10.0])

# β=5.0 confirmation (n=60)

| spec | L=5 | L=7 | L=9 | L=11 |
|---|---|---|---|---|
| flat (1.0, 1.0, 1.0) | 4.58e+08 | 1.57e+09 | 4.58e+09 | 1.22e+10 |
| s4_0.1 (0.1, 1.0, 1.0) | 4.52e+08 | 2.39e+09 | 3.51e+09 | 1.13e+10 |
| s0_0.1 (1.0, 1.0, 0.1) | 1.97e+09 | 1.27e+10 | 2.25e+10 | 4.03e+10 |

ratio s0_0.1/flat: L=5: 4.3, L=7: 8.1, L=9: 4.9, L=11: 3.3
ratio s4_0.1/flat: L=5: 0.99, L=7: 1.52, L=9: 0.77, L=11: 0.93

# H2 classical (OU-modulated diffusion rate), β=4.5, L=9, n=100
| σ | τ_c=0.01 | 1 | 100 | 1e4 | 1e6 | hazard-avg prediction (a=0.6) |
|---|---|---|---|---|---|---|
| 1 | 0.97 | 1.13 | 1.20 | 0.92 | 1.17 | 1.13 |
| 2 | – | 1.60 | 2.12 | – | 2.03 | 1.62 (pure power) |
(values = T_mem / T_mem(flat))

# H2 classical, L-dependence (σ=2, τ_c=100, n=100) vs flat
| L | T_mem (OU) | T_mem (flat) | gain |
|---|---|---|---|
| 5 | 4.48e+07 ± 3.5e+06 | 2.28e+07 | 1.96 |
| 7 | 1.24e+08 ± 1.3e+07 | 7.09e+07 | 1.75 |
| 9 | 2.54e+08 ± 2.7e+07 | 1.20e+08 | 2.12 |
| 11 | 3.41e+08 ± 3.1e+07 | 2.02e+08 | 1.69 |
| 13 | 3.09e+08 ± 2.8e+07 | 2.26e+08 | 1.36 |

# s0 → 0 limit, L=9, n=100 (T/T_flat): s0=0.1: 4.5, 0.01: 29.6, 0.001: 103  — no saturation down to 3 decades

# s0 = 0.01, L-dependence (β=4.5, n=100)
| L | T_mem (s0=0.01) | flat | gain |
|---|---|---|---|
| 5 | 2.97e+08 | 2.28e+07 | 13 |
| 7 | 2.81e+09 | 7.09e+07 | 40 |
| 9 | 3.55e+09 | 1.20e+08 | 30 |
| 11 | 3.92e+09 | 2.02e+08 | 19 |
| 13 | 5.10e+09 | 2.26e+08 | 23 |
Exponent L=5..13: 2.72 (s0=0.01) vs 2.44 (flat); L=7..13 growth: 1.8x (s0=0.01) vs 3.2x (flat). Gain is a prefactor (20-40x), not a change of scaling.

# H2 quantum (qubit gap Δ=4 + one pseudomode, exact Lindblad / Gibbs), scripts/pseudomode_gap.py
- Dynamics (RWA, resonant ω₀=Δ, thermal damping κ, β=2 and 4.5, g=0.05..1, κ=0.03..10): P_exc/P_Boltzmann = 1.0000 for ALL κ, g.
  Non-Markovianity (g>κ) only caps the relaxation rate at κ/2 (golden rule fails), never the population.
- Static mean force (reduced Gibbs state of qubit+mode, dressed-defect population), β=4.5:
  | ω₀ | g=0.1 | g=0.3 | g=0.6 | g=1.0 | 2nd-order Lamb e^{-βg²/(Δ-ω₀)} (g=0.6) |
  |---|---|---|---|---|---|
  | 1 | 0.985 | 0.872 | 0.587 | 0.246 | 0.583 |
  | 2 | 0.978 | 0.820 | 0.473 | 0.155 | 0.445 |
  | 3 | 0.956 | 0.688 | 0.282 | 0.062 | 0.198 |
  | 5 | 1.046 | 1.453 | 3.54 | 16.1 | 5.05 |
  | 6 | 1.023 | 1.219 | 2.11 | 6.45 | 2.25 |
  | 8 | 1.011 | 1.106 | 1.49 | 2.89 | 1.50 |
  Sub-gap bath weight (ω₀<Δ) raises the effective gap and LOWERS the defect population; super-gap weight lowers it.
  Bare-basis 'population' is dominated by virtual admixture of thermally populated sub-gap mode (up to 10^5 × P_B) — not a dynamical defect.
- Artifact documented: counter-rotating coupling + local thermal Lindblad on the mode gives spurious P_exc = κε²/Γ (ε=g/(Δ+ω₀)).

# Lamb-shift lever: flat bath at effective β′ = β(1+δ/Δ), β=4.5, n=100 (scripts/plot_lamb.py, results/lamb_lever.png)
| β′ | δ | exponent (L≤L*) | BH 2.93β′−10.52 | Λ_passive 5→7, 7→9, 9→11, 11→13 | T(13)/T(5) |
|---|---|---|---|---|---|
| 4.5 | 0 | 2.44 | 2.67 | 3.11, 1.69, 1.69, 1.12 | 9.9 |
| 4.70 | 0.18 | 3.35 | 3.25 | 3.06, 2.10, 2.29, 1.64 | 24.0 |
| 4.875 | 0.33 | 3.91 | 3.76 | 3.57, 2.33, 2.77, 1.77 | 40.6 |
| 5.06 | 0.50 | 4.65 | 4.31 | 3.36, 3.19, 2.89, 2.80 | 86.7 |
Λ_passive(11→13) crosses Λ_active = 2.14 (Google 2025) between δ = 0.33 and 0.50 (ω₀ = 2, g ≈ 0.6–1.0).

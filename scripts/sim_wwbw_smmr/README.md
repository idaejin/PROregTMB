# SMMR WW–BW simulation (identification / confounding)

Design: `../../ref/WWBW-SGRQ/SIMULATION-SMMR-DESIGN.md`  
(from vault: `01_PROJECTS/BB/ref/WWBW-SGRQ/SIMULATION-SMMR-DESIGN.md`)

## What is measured

- Bias / RMSE / 95% coverage for between & within slopes  
- Recovery of subject-intercept \(G\) (propensity correlations — **not** visit residual corr)  
- Recovery of \(\phi_\ell\) and \(\tau\) on **unit within-SD** scale  
- Selection by **marginal AIC** (primary) and nℓ (secondary)  
- Roughness proxy for Smooth under Scenario A (apparent curvature)

## Quick

```bash
cd 01_PROJECTS/BB/PROregTMB
Rscript scripts/sim_wwbw_smmr/run_sim.R --quick
```

## Scenario A pilot (recommended next)

```bash
Rscript scripts/sim_wwbw_smmr/run_sim.R --scenario=A --N=200 --nsim=50 --models=M1,M2,M3,M4
```

## Full

```bash
Rscript scripts/sim_wwbw_smmr/run_sim.R --scenario=A --N=500 --nsim=300 --models=M1,M2,M3,M4
Rscript scripts/sim_wwbw_smmr/run_sim.R --scenario=B --N=500 --nsim=300 --models=M1,M2,M3,M4
Rscript scripts/sim_wwbw_smmr/run_sim.R --scenario=C --N=500 --nsim=300 --models=M1,M2,M3,M4,M5
Rscript scripts/sim_wwbw_smmr/run_sim.R --scenario=D --N=500 --nsim=300 --models=M1,M2
```

Outputs: `out/sc*_N*_T*/` (`selection_counts.csv`, `recovery_summary.csv`, `coef_summary.csv`, …).

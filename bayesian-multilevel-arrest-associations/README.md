# Bayesian Models of Recorded Arrest Associations

**How do associations in recorded stop-and-arrest data vary between precincts and times of day?**

This MSc group exercise analyses an administrative dataset using models that allow precincts to differ. It studies conditional associations in the recorded sample; it does not infer a person's criminality or the causal effect of a demographic characteristic.

## How the analysis works

The script constructs first-quarter arrest, demographic, time-of-day, precinct and weapon variables. It fits several Bernoulli logistic models in brms, including varying precinct effects, and compares predictive fit using leave-one-out cross-validation. Posterior predictive checks examine whether the models reproduce features of the observed data.

## Run it

Restore the original `sqf-2025.xlsx` snapshot to this directory. Install readxl, dplyr, brms, loo, writexl and bayesplot, with a working Stan C++ toolchain, then run:

```bash
Rscript analysis.R
```

The workbook and fitted posterior objects are not bundled. Assess MCMC diagnostics and LOO Pareto-k diagnostics before interpreting a new fit. The legacy DIC correction is recorded in [CHANGES.md](../CHANGES.md).

Administrative selection and measurement processes limit interpretation. No causal or policing-fairness conclusion, individual scoring service or independently reproduced numerical result is claimed. The archive does not specify individual ownership; [team authorship](../AUTHORS.md) is retained.

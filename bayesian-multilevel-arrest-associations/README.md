# Bayesian Models of Recorded Arrest Associations

[**Project report**](reports/PROJECT_REPORT.pdf) · [**Review verification**](reports/REVIEW_VERIFICATION.json) · [**Research decision**](docs/RESEARCH_DECISION.md)

This remains a specification/source audit. Original 2025 input, 2024 prior fits, posterior objects and LOO diagnostics are absent. No numerical model ranking or convergence claim is made. An explicit HH:MM/HHMM parser replaces unsafe substring extraction; nine base-R contracts are provided but were not executed in this review.

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

## Time input contract

From this folder, `Rscript parser_checks.R` exercises nine parser examples without fitting a model or requiring brms. It accepts HH:MM or integer HHMM and rejects unsupported encodings. Inspect the workbook before interpreting numeric Excel time fractions. These contracts were not run here because R was unavailable.

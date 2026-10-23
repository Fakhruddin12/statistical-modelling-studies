# Bayesian multilevel arrest associations

Restore the public-source `sqf-2025.xlsx` to this directory, install `readxl`, `dplyr`, `brms`, `loo`, `writexl` and `bayesplot`, and run `Rscript analysis.R` from here. A working Stan C++ toolchain is required. The model uses first-quarter records and constructs arrest, demographic, time-of-day, precinct and weapon variables. The exact workbook snapshot and fitted posterior objects are not bundled.

The code contains multiple Bernoulli-logit models with precinct effects, LOO comparison and posterior predictive checks. Fit diagnostics and Pareto-k diagnostics must be assessed before interpreting a new run. The removed DIC block is explained in the parent change note. No causal or fairness conclusion is inferred from these observational associations.

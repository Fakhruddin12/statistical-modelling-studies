# Reconstructing Infections from Delayed Death Records

[**Project report**](reports/PROJECT_REPORT.pdf) · [**Review verification**](reports/REVIEW_VERIFICATION.json) · [**Research decision**](docs/RESEARCH_DECISION.md)

A recovered 150-day aggregate input supports a new independent Python reconstruction of the archived R equations: 50 smoothing candidates and 200 bootstrap refits. The fatal-infection peak is 20 days before fitted deaths. This is a conditional inverse-problem result, not a forecast of all infections or an R runtime replay.

**If deaths are recorded days or weeks after infection, what can they tell us about when the fatal infections occurred?**

This team project estimates an unobserved infection curve from daily English hospital COVID-19 deaths. It explicitly models the delay between infection and death rather than simply shifting every death backwards by a fixed number of days.

## Method and contribution

A smooth infection curve is represented with B-splines and passed through an assumed infection-to-death delay distribution. A Poisson likelihood compares the resulting expected deaths with the recorded counts. A roughness penalty prevents an excessively wiggly reconstruction; BIC selects its strength. Bootstrap refitting supplies 95% pointwise uncertainty bands.

Fakhruddin Hakim Hussain contributed bootstrap inference, plotting and the report. Mansi Neb and Riccardo Monforte contributed the model matrices, likelihood, gradient and fitting; the source records roughly equal overall effort. [Team record](../AUTHORS.md).

## Run and interpret

Restore `engcov.txt` with columns `nhs` (daily hospital deaths) and `julian` (day index), then run from this directory:

```bash
Rscript analysis.R
```

The code uses base R and splines. The original input is not bundled because redistribution terms are undocumented. Inspect [analysis.R](analysis.R); a fresh numerical reproduction was not performed in this review.

The reconstructed curve concerns **infections leading to the observed hospital deaths**, not all infections in the population. It depends on the delay distribution and smoothing choice. Pointwise bands do not imply simultaneous coverage of the entire curve.

## Independent reference replay

From this folder, install `../review_requirements.txt`, supply the recovered original `engcov.txt`, then run `OPENBLAS_NUM_THREADS=1 python review_replay.py --input engcov.txt --out reports/REVIEW_VERIFICATION.json`. This uses SciPy, not R; optimiser settings and seed are recorded. The reconstructed curve concerns fatal infections linked to these deaths. Private reconstructed curve files are ignored.

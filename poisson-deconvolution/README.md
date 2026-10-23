# Poisson deconvolution

Run `Rscript analysis.R` from this directory after restoring `engcov.txt`. The original input has columns `nhs` (daily hospital deaths) and `julian` (day index). The method convolves a latent B-spline infection trajectory with the specified infection-to-death delay distribution and adds a second-difference penalty. It uses BIC to choose smoothness and bootstrap resampling for uncertainty. Results depend on the assumed delay law and smoothing choice; the inferred curve is not total infections without further ascertainment assumptions.

The supplied data file is omitted from the public package because its redistribution terms are not documented. The script uses base R and `splines`. See the parent README for contributions.

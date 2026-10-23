# Packaging corrections

The supplied Bayesian script estimated DIC by taking column means of pointwise log likelihoods and then summing them. That computes the posterior mean log likelihood, not log likelihood evaluated at posterior-mean parameters. Consequently its effective-parameter calculation collapses to zero up to rounding. The block was removed; the existing LOO comparison is retained. No model was refitted, and no revised numerical results are claimed.

Student numbers were removed from public attribution comments. Original collaborator names and source-recorded contributions are preserved.

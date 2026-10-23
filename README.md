# Statistical Modelling Studies
### Three academic team projects in statistical computing and Bayesian inference

These studies complement the independent research portfolio. They are presented as **collaborative MSc coursework**, with source-level contributions retained.

| Study | Statistical question | Implementation |
|---|---|---|
| [Poisson deconvolution](poisson-deconvolution/) | How can a latent infection curve be reconstructed from delayed death counts? | B-spline basis, delay convolution, second-difference penalty, BIC selection and bootstrap intervals |
| [SEIR contact networks](seir-contact-networks/) | How does household and non-household contact structure change epidemic dynamics? | Stochastic epidemic simulation with explicit contact-network construction |
| [Bayesian multilevel associations](bayesian-multilevel-arrest-associations/) | How do observed arrest associations vary across precincts and time of day? | Bernoulli-logit models in brms, varying precinct effects, posterior checks and LOO comparison |

## Contribution

For the Poisson study, **Fakhruddin Hakim Hussain** contributed bootstrap inference, plotting and the write-up. Mansi Neb contributed model/penalty matrices, likelihood and gradient; Riccardo Monforte contributed fitting and the lambda grid. The source records roughly equal overall effort.

For the SEIR study, the source records approximately equal contributions by the same three authors across implementation, simulation and documentation. The Bayesian group archive does not specify individual ownership, so no sole-authorship claim is made. See [AUTHORS.md](AUTHORS.md).

## Reproduction and evidence

Each folder contains its R script and input notes. Saved numerical results are not available for all studies, so no numerical performance claim is promoted here. This package has not been rerun in R. The Bayesian script uses LOO as its comparison criterion; an incorrect legacy DIC block has been removed and documented in [CHANGES.md](CHANGES.md).

The arrest study is a descriptive conditional-association exercise. The administrative data and observed outcome do not identify causal effects, individual criminality or the fairness of policing decisions. No individual scoring tool is provided.

[Main portfolio](https://github.com/Fakhruddin12)

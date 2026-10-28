# Statistical Modelling Studies

**Three MSc team projects on recovering hidden processes, simulating systems and analysing variation between groups.**

These projects start from different practical questions: when infections happened before recorded deaths, how contact patterns change an epidemic, and how observed arrest associations vary between precincts. Each folder explains the question before the mathematics and code.

| Study | Plain-language question | Technical approach |
|---|---|---|
| [Reconstructing infections from delayed deaths](poisson-deconvolution/) | What infection timeline could have produced the observed deaths? | Poisson deconvolution, penalised splines and bootstrap uncertainty |
| [Epidemics on contact networks](seir-contact-networks/) | How do households and regular contacts affect simulated spread? | Stochastic susceptible–exposed–infectious–removed simulation |
| [Variation in recorded arrest associations](bayesian-multilevel-arrest-associations/) | How do associations vary across precincts and time of day? | Bayesian logistic models with varying precinct effects |

## My contribution and team context

For the infection study, **Fakhruddin Hakim Hussain** contributed bootstrap inference, plotting and the write-up. Mansi Neb contributed model/penalty matrices, likelihood and gradient; Riccardo Monforte contributed fitting and the smoothing grid. The source records roughly equal overall effort.

The SEIR study records approximately equal contributions by the same three authors across implementation, simulation and documentation. The Bayesian archive does not identify individual ownership, so no sole-authorship claim is made. [Contribution record](AUTHORS.md).

## Evidence and reproduction

Each study includes an R script and input/run instructions. Not all original datasets, fitted objects or saved numerical results are available, and these R analyses were not rerun during this documentation review. Numerical success claims are therefore not promoted from missing artifacts.

The Bayesian study uses LOO predictive comparison; the removal of an incorrect legacy DIC block is documented in [CHANGES.md](CHANGES.md). Its observational records support descriptive associations, not causal effects, individual criminality or a verdict on policing fairness.

[Main portfolio](https://github.com/Fakhruddin12)

# Simulating Epidemics with Households and Contact Networks

[**Project report**](reports/PROJECT_REPORT.pdf) · [**Review verification**](reports/REVIEW_VERIFICATION.json) · [**Research decision**](docs/RESEARCH_DECISION.md)

The new Python reference runs 50 replications per case and checks invariants across 200 paths. Its results are new conditional simulations, not restored R figures. Random-only cases raise alpha_random from 0.01 to 0.04, so the design does not isolate network structure.

**Does it matter who meets whom, rather than assuming everyone mixes uniformly?**

This team project simulates infections in a population with households, regular non-household contacts and random contacts. It compares the resulting epidemic paths under different transmission assumptions.

## How it works

People move through susceptible, exposed, infectious and removed states—the SEIR model. The code generates household membership and a contact network, then uses stochastic transitions to simulate spread. Scenarios vary transmission assumptions and whether household/network structure is included.

This is a teaching simulation, not a calibrated public-health forecast. Outcomes depend on the chosen network, transition assumptions and random seed; the repository does not supply a saved benchmark proving one scenario's real-world effect.

## Run it

The script generates its inputs and uses base R. From this directory:

```bash
Rscript analysis.R
```

Set and record a random seed when comparing scenarios. [Source](analysis.R) · [Team contributions](../AUTHORS.md). The source records approximately equal overall contributions by Mansi Neb, Riccardo Monforte and Fakhruddin Hakim Hussain.

## Reproduce the new reference simulations

Install `../review_requirements.txt`, then run `OPENBLAS_NUM_THREADS=1 python review_check.py --replicates 50`. This writes `REVIEW_VERIFICATION.json` in this folder. The published report copy is under `reports/`; new runs should be compared with it explicitly. The original R random stream is not claimed to reproduce.

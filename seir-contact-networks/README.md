# Simulating Epidemics with Households and Contact Networks

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

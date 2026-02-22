# AM_CIF_simulation

## Goal
**Evaluating whether AM-embedded CIFs influence permutation-based estimates of prevalence and heritability.**

To do this, we first generate assortatively mated parental couples and then simulate offspring affection status under two alternative outcome-generating CIFs:

1. **RM-population CIF**: The CIF is calibrated under **randomly mated parents** (`Simulation_CIF_RMpop.r`).
2. **AM-population CIF**: The CIF is calibrated under **assortatively mated parents** (`Simulation_CIF_AMpop.r`).

We then apply the **permutation procedure** described in the paper to break down assortment and compare the resulting **prevalence** and **heritability** estimates across these two scenarios.

## Scenarios
- **Scenario 1 (RM-population CIF)**: `Simulation_CIF_RMpop.r`  
- **Scenario 2 (AM-population CIF)**: `Simulation_CIF_AMpop.r`

In both scenarios, the CIF models offspring risk as a function of:
- **parental affection status**
- **offspring age** (time-to-event / competing-risks framework)

## Output plot
To compare the two scenarios, run:

- `plot Simulation.r`

This script:
1. reads the saved results from both scenarios,
2. combines them into a single dataset, and
3. generates the comparison plot(s).

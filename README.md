# Causal effect estimation from trans-regulatory single-cell CRISPR screens

**Authors:** Oliver P. Christensen, Alex Markham, Hyunseung Kang, Erin Gabriel, and Tune H. Pers.

> **Abstract:** We have reached a pivotal moment in biology: transcriptomic techniques now enable the generation of massive gene expression datasets at single-cell resolution, and CRISPR-based genome editing makes it possible to probe causal relationships in gene regulation. Combining these approaches in single-cell CRISPR screens allows researchers to move beyond correlation and begin to reason about causality. Going forward, such causal relationships between perturbations and genes are likely to greatly deepen our understanding of cellular functions and disease mechanisms. However, estimating causal effects requires additional causal assumptions, which are needed to give estimates from traditional statistical methods a causal interpretation. Here, we first introduce key concepts and principles of causal effect estimation in the analysis of trans-regulatory single-cell CRISPRi data. We then outline a set of causal assumptions, which allow for a causal interpretation of estimates from existing statistical methods that can be used to estimate associations between parturbations and genes from single-cell CRISPR screen data. Additionally, we also briefly give an overview of the statistical methods. Finally, we provide an illustrative example to showcase how the causal effects can be biased when any of these assumptions are violated.

---

# Worked Example: Violations of Causal Assumptions in a Simple Simulation

## Overview

This simulation study illustrates how violations of key causal inference assumptions can bias estimated treatment effects in a setting resembling a trans-regulatory single-guide single-cell CRISPRi screen.

We consider:

- A binary perturbation indicator `t` (gRNA targeting gene X present vs absent)
- A count outcome `Y` (expression of gene Y)
- A Poisson data-generating process
- True causal log fold change = −1 (fold change ≈ 0.368)

For each scenario, data are simulated, a Poisson regression model (`Y ~ t`) is fitted, and estimates are averaged across Monte Carlo replicates.

---

## Scenarios

### Baseline (No Assumption Violations)

Treatment assignment is randomized, and there are no confounders, effect modifiers, or heterogeneous treatment versions.

All causal assumptions are satisfied, and the estimated effect recovers the true causal fold change.

---

### Consistency Violation

Multiple gRNAs targeting gene X have different knockdown efficiencies but are combined into a single treatment category.

Specifically:

- gRNA1: moderate effect
- gRNA2: strong effect
- gNT1: non-targeting control

Aggregating these guides into one indicator means that “treatment” does not correspond to a single well-defined intervention, leading to biased estimates.

---

### No Interference Violation

An environmental variable modifies the treatment effect.

The outcome depends on the interaction between treatment and environment:

    effect = beta1 × t × env

This creates heterogeneous treatment effects across cells, violating the assumption that outcomes depend only on a unit’s own treatment.

---

### Conditional Ignorability Violation From Confounder

An unobserved binary variable C:

- Increases the probability of receiving treatment
- Directly increases the outcome Y

Because C is omitted from the analysis model, treatment assignment is confounded, producing biased effect estimates.

---

### Conditional Ignorability Violation From Collider

An observed variable C:
- Is affected by both treatment and outcome Y

Because the analysis conditions on C, this opens a backdoor path and induces bias in the estimated treatment effect.

---

## Output

For each scenario, the script reports:

- Estimated baseline mean expression
- Estimated treatment fold change

---

## Requirements

R packages:

- tidyverse
- broom

Install if needed:

    install.packages(c("tidyverse", "broom"))

---

## Reproducibility

A fixed random seed is used to ensure reproducible results.

---

## Interpretation

The simulations demonstrate that even simple violations of causal assumptions can produce substantial bias in estimated effects, emphasizing the importance of carefully considering these assumptions when interpreting results from CRISPR screening experiments or similar observational analyses.

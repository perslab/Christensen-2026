############################################################
# Simulation study illustrating violations of causal assumptions
# in a simple Poisson regression setting.
#
# Context:
# This script mimics a low-MOI single-cell CRISPRi screen
# with a binary perturbation indicator (t) for gene X and a
# count outcome (Y) representing expression of gene Y.
#
# Across all scenarios, the true causal effect of perturbing X
# on Y is a log fold change of -1 (fold change ≈ 0.368).
#
# For each scenario:
#   1) Data are simulated under a specified mechanism
#   2) A Poisson regression Y ~ t is fitted
#   3) Estimated coefficients are stored
#   4) Results are averaged across Monte Carlo replicates
#
# The scenarios correspond to:
#   - No assumption violations
#   - Violation of conditional ignorability
#   - Violation of no interference
#   - Violation of consistency
############################################################

library(tidyverse)
library(broom)

set.seed(1)  # Ensures reproducibility


############################################################
# 1) BASELINE: NO ASSUMPTION VIOLATIONS
#
# Treatment is randomized and there are no confounders,
# effect modifiers, or multiple treatment versions.
############################################################

estimates_no_violation <- data.frame(beta0 = NA, beta1 = NA)

for (i in 1:100) {
  
  n <- 10000            # Number of cells
  p <- 0.5              # Probability of receiving treatment
  t <- rbinom(n, 1, p)  # Randomized treatment assignment
  
  beta0 <- 2            # Baseline log mean expression
  beta1 <- -1           # True log fold change (treatment effect)
  
  # Outcome generated from Poisson model
  Y <- rpois(n, lambda = exp(beta0 + beta1 * t))
  
  df <- data.frame(Y = Y, t = t)
  
  # Estimate effect using Poisson regression
  pois_model <- glm(Y ~ t, family = "poisson", data = df)
  tidy_sum <- tidy(pois_model)
  
  estimates_no_violation[i, ] <- c(tidy_sum$estimate[1], tidy_sum$estimate[2])
}

############################################################
# 2) CONSISTENCY VIOLATION
#
# Multiple gRNAs targeting gene X have different knockdown
# efficiencies but are aggregated into a single treatment
# indicator. Thus, "treatment" does not represent a single
# well-defined intervention.
############################################################

estimates_cons <- data.frame(beta0 = NA, beta1 = NA)

for (i in 1:100) {
  
  n <- 10000
  
  # Three guide types:
  #   gNT1  : non-targeting control
  #   gRNA1 : moderate knockdown
  #   gRNA2 : strong knockdown
  t <- sample(c("gRNA1", "gRNA2", "gNT1"), n, TRUE)
  
  # Log effects for each guide
  lambda_map <- c(
    gNT1  = 0,
    gRNA1 = -0.5,
    gRNA2 = -1
  )
  
  beta0 <- 2
  
  # Outcome depends on specific guide used
  Y <- rpois(n, lambda = exp(beta0 + lambda_map[t]))
  
  # Aggregate targeting guides into one category
  t_aggregate <- t
  t_aggregate[grepl("gRNA", t)] <- "gRNA"
  
  df <- data.frame(Y = Y, t = t_aggregate)
  
  pois_model <- glm(Y ~ t, family = "poisson", data = df)
  tidy_sum <- tidy(pois_model)
  
  estimates_cons[i, ] <- c(tidy_sum$estimate[1], tidy_sum$estimate[2])
}

############################################################
# 3) NO INTERFERENCE VIOLATION
#
# An environmental variable modifies the treatment effect.
# This produces heterogeneous treatment effects across cells.
############################################################

estimates_no_int <- data.frame(beta0 = NA, beta1 = NA)

for (i in 1:100) {
  
  n <- 10000
  p <- 0.5
  t <- rbinom(n, 1, p)
  
  # Environmental factor affecting treatment response
  env <- rgamma(n, shape = 0.9)
  
  beta0 <- 2
  beta1 <- -1
  
  # Treatment effect varies with environment
  Y <- rpois(n, exp(beta0 + beta1 * t * env))
  
  df <- data.frame(Y = Y, t = t)
  
  pois_model <- glm(Y ~ t, family = "poisson", data = df)
  tidy_sum <- tidy(pois_model)
  
  estimates_no_int[i, ] <- c(tidy_sum$estimate[1], tidy_sum$estimate[2])
}

############################################################
# 4) CONDITIONAL IGNORABILITY VIOLATION FROM CONFOUNDER
#
# An unobserved confounder C affects both:
#   - Probability of treatment
#   - Outcome Y
#
# The analysis model omits C, inducing bias.
############################################################

estimates_cond_ignor_con <- data.frame(beta0 = NA, beta1 = NA)

for (i in 1:100) {
  
  n <- 10000
  
  pC <- 0.4
  C <- rbinom(n, 1, pC)        # Unobserved confounder
  
  # Treatment probability depends on C
  p <- 0.5 + 0.3 * C
  t <- rbinom(n, 1, p)
  
  beta0 <- 2
  beta1 <- -1                  # True causal effect
  beta2 <- 1                   # Effect of confounder on outcome
  
  # Outcome depends on both treatment and confounder
  Y <- rpois(n, lambda = exp(beta0 + beta1 * t + beta2 * C))
  
  df <- data.frame(Y = Y, t = t, C = C)
  
  # Confounder omitted from model
  pois_model <- glm(Y ~ t, family = "poisson", data = df)
  tidy_sum <- tidy(pois_model)
  
  estimates_cond_ignor_con[i, ] <- c(tidy_sum$estimate[1], tidy_sum$estimate[2])
}


############################################################
# 5) CONDITIONAL IGNORABILITY VIOLATION FROM COLLIDER
#
# An observed collider C is affected by both:
#   - Probability of treatment
#   - Outcome Y
#
# The analysis model condition on C, inducing bias.
############################################################

estimates_cond_ignor_col <- data.frame(beta0 = NA, beta1 = NA)

for (i in 1:100) {
  
  n <- 10000            # Number of cells
  p <- 0.5              # Probability of receiving treatment
  t <- rbinom(n, 1, p)  # Randomized treatment assignment
  
  beta0 <- 2            # Baseline log mean expression
  beta1 <- -1           # True log fold change (treatment effect)
  beta2 <- 3            # Effect of treatment on collider
  beta3 <- -2            # Effect of outcome on collider
  
  # Outcome generated from Poisson model
  Y <- rpois(n, lambda = exp(beta0 + beta1 * t))
  
  # Collider generated from Poisson model as a function of both treatment and outcome
  C <- rpois(n, lambda = exp(beta0 + beta2 * t + beta3*Y))
  
  df <- data.frame(Y = Y, t = t, C = C)
  
  # Collider included in the model
  pois_model <- glm(Y ~ t + C, family = "poisson", data = df)
  tidy_sum <- tidy(pois_model)
  
  estimates_cond_ignor_col[i, ] <- c(tidy_sum$estimate[1], tidy_sum$estimate[2])
}


# Calculate statistics for each scenario

# Helper function to compute summary statistics
summarize_estimates <- function(est_df) {
  c(
    mean_beta0 = mean(exp(est_df$beta0)),
    mean_beta1 = mean(exp(est_df$beta1))
  )
}

# Combine into a data.frame
results <- data.frame(
  "No Violation" = summarize_estimates(estimates_no_violation),
  "Consistency" = summarize_estimates(estimates_cons),
  "No Interference"  = summarize_estimates(estimates_no_int),
  "Conditional Ignorability Confounder" = summarize_estimates(estimates_cond_ignor_con),
  "Conditional Ignorability Collider" = summarize_estimates(estimates_cond_ignor_col)
)

rownames(results) <- c("Mean beta0", "Mean beta1")

# View results
print(results)

# Save as CSV
write.csv(results, "causal_simulation_results.csv", row.names = FALSE)
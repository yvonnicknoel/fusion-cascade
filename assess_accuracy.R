## Assess true/false positive rates for Table 7 scenarios
## This shows how many TRUE differences each method detects

library(parallel)
library(multcomp)

source("power_simulations.R")

# Set parameters
n <- 20
sigma <- 1
alpha <- 0.05
R <- 2000

set.seed(123)  # For reproducibility

cat("Assessing accuracy (TP/FP) for Table 7 scenarios (n=20)\n")
cat("R =", R, "replications per scenario\n\n")

## Scenario 1: Monotonic k=5, delta=0.5
cat("=" , rep("=", 60), "\n", sep="")
cat("Scenario: Monotonic k=5, delta=0.5\n")
cat("True structure: All 5 groups different -> 10 true pairwise differences\n")
mu1 <- scenario_monotonic(k = 5, effect_size = 0.5, sigma = sigma)
c1 <- compute_c_approx(k=5)
power1 <- assess_power(mu1, n = n, sigma = sigma, c_factor = c1,
                       alpha = alpha, R = R, compute_accuracy = TRUE)
cat("\nTrue differences:", power1$n_true_diff[1], "\n")
cat("True non-differences:", power1$n_true_same[1], "\n\n")
print(power1[, c("method", "mean_rejections", "mean_tp", "mean_fp")])
cat("\n")

## Scenario 2: Monotonic k=5, delta=1.0
cat("=" , rep("=", 60), "\n", sep="")
cat("Scenario: Monotonic k=5, delta=1.0\n")
cat("True structure: All 5 groups different -> 10 true pairwise differences\n")
mu2 <- scenario_monotonic(k = 5, effect_size = 1.0, sigma = sigma)
c2 <- compute_c_approx(k=5)
power2 <- assess_power(mu2, n = n, sigma = sigma, c_factor = c2,
                       alpha = alpha, R = R, compute_accuracy = TRUE)
cat("\nTrue differences:", power2$n_true_diff[1], "\n")
cat("True non-differences:", power2$n_true_same[1], "\n\n")
print(power2[, c("method", "mean_rejections", "mean_tp", "mean_fp")])
cat("\n")

## Scenario 3: Monotonic k=5, delta=2.0
cat("=" , rep("=", 60), "\n", sep="")
cat("Scenario: Monotonic k=5, delta=2.0\n")
cat("True structure: All 5 groups different -> 10 true pairwise differences\n")
mu3 <- scenario_monotonic(k = 5, effect_size = 2.0, sigma = sigma)
c3 <- compute_c_approx(k=5)
power3 <- assess_power(mu3, n = n, sigma = sigma, c_factor = c3,
                       alpha = alpha, R = R, compute_accuracy = TRUE)
cat("\nTrue differences:", power3$n_true_diff[1], "\n")
cat("True non-differences:", power3$n_true_same[1], "\n\n")
print(power3[, c("method", "mean_rejections", "mean_tp", "mean_fp")])
cat("\n")

## Scenario 4: TwoGroups k=6, delta=0.5
cat("=" , rep("=", 60), "\n", sep="")
cat("Scenario: TwoGroups k=6 (3+3), delta=0.5\n")
cat("True structure: {1,2,3} vs {4,5,6} -> 9 true pairwise differences\n")
mu4 <- scenario_two_groups(k = 6, k1 = 3, effect_size = 0.5, sigma = sigma)
c4 <- compute_c_approx(k=6)
power4 <- assess_power(mu4, n = n, sigma = sigma, c_factor = c4,
                       alpha = alpha, R = R, compute_accuracy = TRUE)
cat("\nTrue differences:", power4$n_true_diff[1], "\n")
cat("True non-differences:", power4$n_true_same[1], "\n\n")
print(power4[, c("method", "mean_rejections", "mean_tp", "mean_fp")])
cat("\n")

## Scenario 5: TwoGroups k=6, delta=1.0
cat("=" , rep("=", 60), "\n", sep="")
cat("Scenario: TwoGroups k=6 (3+3), delta=1.0\n")
cat("True structure: {1,2,3} vs {4,5,6} -> 9 true pairwise differences\n")
mu5 <- scenario_two_groups(k = 6, k1 = 3, effect_size = 1.0, sigma = sigma)
c5 <- compute_c_approx(k=6)
power5 <- assess_power(mu5, n = n, sigma = sigma, c_factor = c5,
                       alpha = alpha, R = R, compute_accuracy = TRUE)
cat("\nTrue differences:", power5$n_true_diff[1], "\n")
cat("True non-differences:", power5$n_true_same[1], "\n\n")
print(power5[, c("method", "mean_rejections", "mean_tp", "mean_fp")])
cat("\n")

## Scenario 6: TwoGroups k=6, delta=2.0
cat("=" , rep("=", 60), "\n", sep="")
cat("Scenario: TwoGroups k=6 (3+3), delta=2.0\n")
cat("True structure: {1,2,3} vs {4,5,6} -> 9 true pairwise differences\n")
mu6 <- scenario_two_groups(k = 6, k1 = 3, effect_size = 2.0, sigma = sigma)
c6 <- compute_c_approx(k=6)
power6 <- assess_power(mu6, n = n, sigma = sigma, c_factor = c6,
                       alpha = alpha, R = R, compute_accuracy = TRUE)
cat("\nTrue differences:", power6$n_true_diff[1], "\n")
cat("True non-differences:", power6$n_true_same[1], "\n\n")
print(power6[, c("method", "mean_rejections", "mean_tp", "mean_fp")])
cat("\n")

# Combine results for export
all_results <- rbind(
  cbind(scenario = "Monotonic k=5 delta=0.5", power1[, c("method", "mean_rejections", "mean_tp", "mean_fp")]),
  cbind(scenario = "Monotonic k=5 delta=1.0", power2[, c("method", "mean_rejections", "mean_tp", "mean_fp")]),
  cbind(scenario = "Monotonic k=5 delta=2.0", power3[, c("method", "mean_rejections", "mean_tp", "mean_fp")]),
  cbind(scenario = "TwoGroups k=6 delta=0.5", power4[, c("method", "mean_rejections", "mean_tp", "mean_fp")]),
  cbind(scenario = "TwoGroups k=6 delta=1.0", power5[, c("method", "mean_rejections", "mean_tp", "mean_fp")]),
  cbind(scenario = "TwoGroups k=6 delta=2.0", power6[, c("method", "mean_rejections", "mean_tp", "mean_fp")])
)

write.csv(all_results, "accuracy_assessment.csv", row.names = FALSE)
cat("\n=== Results saved to accuracy_assessment.csv ===\n")

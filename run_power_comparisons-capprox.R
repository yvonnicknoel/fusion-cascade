## ------------------------------------------------------------
## Run complete power comparison study
## Uses calibrated c values from the CSV file
## ------------------------------------------------------------

library(parallel)

source("power_simulations.R")

## ------------------------------------------------------------
## 1. Define scenarios to test
## ------------------------------------------------------------

scenarios_list <- list(
  list(
    name = "Monotonic_k5_d0.5",
    k = 5,
    mu_fn = function() scenario_monotonic(k = 5, effect_size = 0.5)
  ),
  list(
    name = "Monotonic_k5_d1.0",
    k = 5,
    mu_fn = function() scenario_monotonic(k = 5, effect_size = 1.0)
  ),
  list(
    name = "TwoGroups_k6_d1.0",
    k = 6,
    mu_fn = function() scenario_two_groups(k = 6, k1 = 3, effect_size = 1.0)
  ),
  list(
    name = "TwoGroups_k6_d0.5",
    k = 6,
    mu_fn = function() scenario_two_groups(k = 6, k1 = 3, effect_size = 0.5)
  ),
  list(
    name = "ThreeGroups_k9",
    k = 9,
    mu_fn = function() scenario_three_groups(k = 9, k1 = 3, k2 = 3,
                                              effect_size1 = 0.5,
                                              effect_size2 = 1.0)
  ),
  list(
    name = "Outlier_k5_d2.0",
    k = 5,
    mu_fn = function() scenario_outlier(k = 5, effect_size = 2.0)
  ),
  list(
    name = "Outlier_k5_d1.0",
    k = 5,
    mu_fn = function() scenario_outlier(k = 5, effect_size = 1.0)
  ),
  list(
    name = "TwoGroups_k8_d0.8",
    k = 8,
    mu_fn = function() scenario_two_groups(k = 8, k1 = 4, effect_size = 0.8)
  )
)

## ------------------------------------------------------------
## 2. Run power assessment for one (scenario, n) combination
## ------------------------------------------------------------

run_one_power <- function(scenario, n, sigma = 1, alpha = 0.05, R = 2000) {
  k <- scenario$k
  scenario_name <- scenario$name
  mu <- scenario$mu_fn()

  # Get calibrated c for this (k, n)
  # calib_row <- calib_data[calib_data$k == k & calib_data$n == n, ]

  # if (nrow(calib_row) == 0) {
  #   cat("Warning: No calibration for k =", k, "n =", n, "\n")
  #   return(NULL)
  # }

  # c_factor <- calib_row$c[1]

  # Compute c factor following formula (9)
  c_factor = compute_c_approx(k)

  cat("Running:", scenario_name, "k =", k, "n =", n, "c =", round(c_factor, 6), "\n")

  power_res <- assess_power(mu, n = n, sigma = sigma, c_factor = c_factor,
                            alpha = alpha, R = R)

  power_res$scenario <- scenario_name
  power_res$k <- k
  power_res$n <- n
  power_res$c_factor <- c_factor

  power_res
}

## ------------------------------------------------------------
## 3. Grid of (scenario, n) combinations
## ------------------------------------------------------------

# Sample sizes to test
ns_to_test <- c(10, 15, 20, 25, 30)

# Build grid
power_grid <- expand.grid(
  scenario_idx = seq_along(scenarios_list),
  n = ns_to_test,
  stringsAsFactors = FALSE
)

## ------------------------------------------------------------
## 4. Run in parallel
## ------------------------------------------------------------

n_cores <- max(1L, detectCores() - 2L)
set.seed(456)

run_one_row <- function(idx) {
  scenario <- scenarios_list[[power_grid$scenario_idx[idx]]]
  n <- power_grid$n[idx]

  run_one_power(scenario, n = n, sigma = 1, alpha = 0.05, R = 2000)
}

power_results_list <- mclapply(
  X = seq_len(nrow(power_grid)),
  FUN = run_one_row,
  mc.cores = n_cores
)

# Remove NULLs
power_results_list <- power_results_list[!sapply(power_results_list, is.null)]

# Combine
power_results <- do.call(rbind, power_results_list)

# Reorder columns
power_results <- power_results[, c("scenario", "k", "n", "c_factor", "method",
                                   "power", "mean_rejections")]

## Save results
write.csv(power_results, "power_comparison_results_capprox.csv", row.names = FALSE)

cat("\n=== Power comparison complete ===\n")
cat("Results saved to power_comparison_results.csv\n")

## ------------------------------------------------------------
## 5. Quick summary
## ------------------------------------------------------------

cat("\n=== Summary by scenario (averaged over n) ===\n")
library(dplyr)

# Set method order: Fusion, Tukey, BH, Holm, Scheffe
method_order <- c("Fusion", "Tukey", "BH", "Holm", "Scheffe")

summary_by_scenario <- power_results %>%
  group_by(scenario, method) %>%
  summarise(
    mean_power = mean(power),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(names_from = method, values_from = mean_power) %>%
  select(scenario, all_of(method_order))

print(summary_by_scenario)

cat("\n=== Summary by method (averaged over all scenarios) ===\n")
summary_by_method <- power_results %>%
  group_by(method) %>%
  summarise(
    mean_power = mean(power),
    sd_power = sd(power),
    .groups = "drop"
  ) %>%
  mutate(method = factor(method, levels = method_order)) %>%
  arrange(method)

print(summary_by_method)

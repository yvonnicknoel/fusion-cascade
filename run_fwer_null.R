## FWER assessment under the global null
## All methods compared with all group means equal

source("power_simulations.R")

cat("Running FWER simulations under global null...\n")

# Load calibration data to get c values
# calib_data <- read.csv("nested_models_calibration_fusion.csv")

# Grid of (k, n) combinations to test
k_values <- c(5, 6, 8, 10, 12, 15)
n_values <- c(10, 15, 20, 25, 30)

# Create grid
fwer_grid <- expand.grid(k = k_values, n = n_values)

# Number of replications
R <- 2000
alpha <- 0.05

results_list <- list()

for (idx in 1:nrow(fwer_grid)) {
  k <- fwer_grid$k[idx]
  n <- fwer_grid$n[idx]

  cat("Testing k =", k, ", n =", n, "\n")

  # Get calibrated c for fusion method
  # c_row <- calib_data[calib_data$k == k & calib_data$n == n, ]
  # if (nrow(c_row) == 0) {
  #   cat("  Warning: No calibration for k =", k, ", n =", n, "; skipping\n")
  #   next
  # }
  # c_factor <- c_row$c[1]

  # Use the generic approximate formula for c
  c_factor = compute_c_approx(k)

  # Global null: all means equal to zero
  mu_null <- rep(0, k)

  # Run power assessment (which under null gives us FWER)
  set.seed(1000 + idx)
  result <- assess_power(mu_null, n = n, sigma = 1,
                        c_factor = c_factor, alpha = alpha, R = R)

  # Add k and n to results
  result$k <- k
  result$n <- n
  result$c_factor <- c_factor

  results_list[[idx]] <- result
}

# Combine results
fwer_results <- do.call(rbind, results_list)

# Reorder columns
fwer_results <- fwer_results[, c("k", "n", "method", "power", "c_factor")]
names(fwer_results)[names(fwer_results) == "power"] <- "fwer"

# Save results
write.csv(fwer_results, "fwer_null_results.csv", row.names = FALSE)

cat("\n=== FWER under global null completed ===\n")
cat("Results saved to fwer_null_results.csv\n\n")

# Print summary table for selected (k, n) combinations
cat("=== Summary for manuscript (selected combinations) ===\n\n")

selected_combos <- data.frame(
  k = c(5, 5, 8, 8, 10, 10),
  n = c(10, 20, 10, 20, 10, 20)
)

for (i in 1:nrow(selected_combos)) {
  k_sel <- selected_combos$k[i]
  n_sel <- selected_combos$n[i]

  subset_data <- fwer_results[fwer_results$k == k_sel & fwer_results$n == n_sel, ]

  if (nrow(subset_data) > 0) {
    cat(sprintf("k=%d, n=%d:\n", k_sel, n_sel))
    for (j in 1:nrow(subset_data)) {
      cat(sprintf("  %s: %.4f\n",
                  subset_data$method[j],
                  subset_data$fwer[j]))
    }
    cat("\n")
  }
}

# Overall summary by method
cat("\n=== Overall FWER by method (averaged across all k,n) ===\n")
library(dplyr)
overall <- fwer_results %>%
  group_by(method) %>%
  summarise(
    mean_fwer = mean(fwer),
    sd_fwer = sd(fwer),
    min_fwer = min(fwer),
    max_fwer = max(fwer),
    .groups = "drop"
  )
print(overall)

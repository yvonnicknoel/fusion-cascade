## Interpretational Coherence Analysis
## Compare partition recovery and transitivity between Fusion and Tukey HSD
## Tests two scenarios: moderate and strong effect sizes

source("power_simulations.R")

cat("Running interpretational coherence analysis...\n")

# Common parameters
k <- 5
n <- 20
sigma <- 1
alpha <- 0.05
R <- 1000

# Two scenarios
scenarios <- list(
  moderate = list(
    mu = c(0, 0, 0.5, 1.0, 1.0),
    label = "Moderate effects"
  ),
  strong = list(
    mu = c(0, 0, 1.0, 2.0, 2.0),
    label = "Strong effects"
  )
)

cat("Common parameters: k =", k, ", n =", n, ", R =", R, "\n")
cat("True partition for both: {{1,2}, {3}, {4,5}}\n\n")

# Get calibrated c for fusion method
c_factor <- compute_c_approx(k)
cat("Using calibration factor c =", c_factor, "\n\n")

# Helper function to convert partition to string for comparison
partition_to_string <- function(partition) {
  # Sort each group, then sort groups
  sorted_partition <- lapply(partition, sort)
  sorted_partition <- sorted_partition[order(sapply(sorted_partition, function(x) x[1]))]
  paste(sapply(sorted_partition, function(grp) paste("{", paste(grp, collapse=","), "}", sep="")),
        collapse="")
}

# Define the true partition as string
true_partition_str <- partition_to_string(list(c(1,2), c(3), c(4,5)))

# Storage for results across scenarios
all_results <- list()

# Helper function to check if Tukey pattern is transitive
is_transitive <- function(rejection_matrix) {
  # rejection_matrix[i,j] = TRUE if we reject H0: mu_i = mu_j
  k <- nrow(rejection_matrix)

  # Check transitivity: if i~j and j~k, then i~k
  # (where ~ means "not significantly different")
  for (i in 1:(k-2)) {
    for (j in (i+1):(k-1)) {
      for (m in (j+1):k) {
        # If i~j (not rejected) and j~m (not rejected), then i~m should not be rejected
        if (!rejection_matrix[i,j] && !rejection_matrix[j,m] && rejection_matrix[i,m]) {
          return(FALSE)
        }
        # If i!=j (rejected) and j!=m (rejected), then i!=m should be rejected
        # This is the other direction of transitivity
        if (rejection_matrix[i,j] && rejection_matrix[j,m] && !rejection_matrix[i,m]) {
          return(FALSE)
        }
      }
    }
  }
  return(TRUE)
}

set.seed(2000)

# Loop over scenarios
for (scenario_name in names(scenarios)) {
  scenario <- scenarios[[scenario_name]]
  mu_true <- scenario$mu

  cat("==============================================\n")
  cat("Scenario:", scenario$label, "\n")
  cat("True means:", mu_true, "\n")
  cat("==============================================\n\n")

  # Storage for results
  fusion_partitions <- vector("list", R)
  fusion_exact_match <- 0
  fusion_merged_123 <- 0
  fusion_other <- 0

  tukey_consistent <- 0
  tukey_nontransitive <- 0

for (r in 1:R) {
  if (r %% 100 == 0) cat("Replication", r, "/", R, "\n")

  # Simulate data
  dat <- simulate_data(mu_true, n, sigma)

  # --- FUSION CASCADE ---
  # Call cascade_fusions directly to get detailed test results
  m <- k - 1
  alpha_seq <- c_factor * alpha / (m:1)
  fusion_result <- cascade_fusions(dat, alpha_seq = alpha_seq)

  # Reconstruct partition from fusion tests
  # Start with each group in its own cluster
  partition <- as.list(1:k)

  if (!is.null(fusion_result$tests)) {
    # Apply fusions sequentially until rejection
    for (i in 1:nrow(fusion_result$tests)) {
      test <- fusion_result$tests[i, ]

      # If this test was not rejected, apply the fusion
      if (test$p >= test$alpha) {
        # Find which clusters contain g1 and g2
        g1 <- as.numeric(test$g1)
        g2 <- as.numeric(test$g2)

        idx1 <- which(sapply(partition, function(x) g1 %in% x))
        idx2 <- which(sapply(partition, function(x) g2 %in% x))

        if (length(idx1) > 0 && length(idx2) > 0 && idx1 != idx2) {
          # Merge the two clusters
          partition[[idx1]] <- c(partition[[idx1]], partition[[idx2]])
          partition[[idx2]] <- NULL
        }
      } else {
        # Rejection: stop here
        break
      }
    }
  }

  fusion_partitions[[r]] <- partition
  partition_str <- partition_to_string(partition)

  # Check if it matches true partition exactly
  if (partition_str == true_partition_str) {
    fusion_exact_match <- fusion_exact_match + 1
  } else if (partition_str == partition_to_string(list(c(1,2,3), c(4,5)))) {
    # Check if groups 1,2,3 merged together
    fusion_merged_123 <- fusion_merged_123 + 1
  } else {
    fusion_other <- fusion_other + 1
  }

  # --- TUKEY HSD ---
  tukey_result <- tukey_posthoc(dat, alpha = alpha)

  # Build rejection matrix
  rejection_matrix <- matrix(FALSE, k, k)
  if (tukey_result$n_rejections > 0) {
    # Parse the names from TukeyHSD output (format: "2-1", "3-1", etc.)
    rejection_names <- names(tukey_result$rejections)[tukey_result$rejections]
    for (name in rejection_names) {
      parts <- as.numeric(strsplit(name, "-")[[1]])
      g1 <- parts[1]
      g2 <- parts[2]
      rejection_matrix[g1, g2] <- TRUE
      rejection_matrix[g2, g1] <- TRUE
    }
  }

  # Check expected pattern for true partition {{1,2}, {3}, {4,5}}
  # Should reject: 1-3, 2-3, 3-4, 3-5
  # Should not reject: 1-2, 4-5
  expected_pattern <- rejection_matrix[1,3] && rejection_matrix[2,3] &&
                      rejection_matrix[3,4] && rejection_matrix[3,5] &&
                      !rejection_matrix[1,2] && !rejection_matrix[4,5]

  if (expected_pattern) {
    tukey_consistent <- tukey_consistent + 1
  }

  # Check transitivity
  if (!is_transitive(rejection_matrix)) {
    tukey_nontransitive <- tukey_nontransitive + 1
  }
}

  # Print results for this scenario
  cat("\n--- Results for", scenario$label, "---\n\n")

  cat("FUSION CASCADE:\n")
  cat(sprintf("  Exact true partition {{1,2}, {3}, {4,5}}: %.1f%%\n",
              100 * fusion_exact_match / R))
  cat(sprintf("  Merged to {{1,2,3}, {4,5}}: %.1f%%\n",
              100 * fusion_merged_123 / R))
  cat(sprintf("  Other partitions: %.1f%%\n",
              100 * fusion_other / R))
  cat("\n")

  cat("TUKEY HSD:\n")
  cat(sprintf("  Consistent with true partition: %.1f%%\n",
              100 * tukey_consistent / R))
  cat(sprintf("  Non-transitive patterns: %.1f%%\n",
              100 * tukey_nontransitive / R))
  cat("\n\n")

  # Store results
  all_results[[scenario_name]] <- data.frame(
    scenario = scenario$label,
    effect_size = ifelse(scenario_name == "moderate", "$\\delta = 0.5$", "$\\delta = 1.0$"),
    fusion_exact = 100 * fusion_exact_match / R,
    fusion_other = 100 * fusion_other / R,
    tukey_consistent = 100 * tukey_consistent / R,
    tukey_nontransitive = 100 * tukey_nontransitive / R
  )
}

# Combine all results
final_results <- do.call(rbind, all_results)
rownames(final_results) <- NULL

cat("\n=== SUMMARY TABLE ===\n\n")
print(final_results)

# Save results
write.csv(final_results, "interpretational_coherence_results_capprox.csv", row.names = FALSE)
cat("\nResults saved to interpretational_coherence_results_capprox.csv\n")

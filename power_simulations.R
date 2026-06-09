## ------------------------------------------------------------
## Power simulations: Comparison of fusion cascade vs classical methods
## ------------------------------------------------------------
## Compares:
##  - Calibrated fusion cascade
##  - Tukey HSD
##  - Holm-Bonferroni (step-down)
##  - Benjamini-Hochberg (FDR control)
##  - Scheffé method
## ------------------------------------------------------------

library(parallel)
library(multcomp)  # for Tukey HSD

source("nested_models.R")

## ------------------------------------------------------------
## 1. Classical post-hoc methods
## ------------------------------------------------------------

## Tukey HSD
tukey_posthoc <- function(dat, alpha = 0.05) {
  fit <- aov(y ~ group, data = dat)
  tukey_result <- TukeyHSD(fit, conf.level = 1 - alpha)

  # Extract p-values for all pairwise comparisons
  p_values <- tukey_result$group[, "p adj"]
  rejections <- p_values < alpha

  list(
    any_reject = any(rejections),
    n_rejections = sum(rejections),
    p_values = p_values,
    rejections = rejections
  )
}

## Holm-Bonferroni on all pairwise comparisons
holm_posthoc <- function(dat, alpha = 0.05) {
  k <- nlevels(dat$group)
  means <- tapply(dat$y, dat$group, mean)
  lev <- names(means)

  # All pairwise t-tests
  pairs <- t(combn(lev, 2))
  n_pairs <- nrow(pairs)

  # Pooled error from full ANOVA
  fit <- lm(y ~ group, data = dat)
  mse <- sum(residuals(fit)^2) / df.residual(fit)
  n_per_group <- table(dat$group)

  p_values <- numeric(n_pairs)
  for (i in 1:n_pairs) {
    g1 <- pairs[i, 1]
    g2 <- pairs[i, 2]

    diff <- abs(means[g1] - means[g2])
    se <- sqrt(mse * (1/n_per_group[g1] + 1/n_per_group[g2]))
    t_stat <- diff / se
    p_values[i] <- 2 * pt(-abs(t_stat), df = df.residual(fit))
  }

  # Holm adjustment
  p_adjusted <- p.adjust(p_values, method = "holm")
  rejections <- p_adjusted < alpha

  list(
    any_reject = any(rejections),
    n_rejections = sum(rejections),
    p_values = p_adjusted,
    rejections = rejections
  )
}

## Benjamini-Hochberg (FDR control)
bh_posthoc <- function(dat, alpha = 0.05) {
  k <- nlevels(dat$group)
  means <- tapply(dat$y, dat$group, mean)
  lev <- names(means)

  # All pairwise t-tests
  pairs <- t(combn(lev, 2))
  n_pairs <- nrow(pairs)

  # Pooled error from full ANOVA
  fit <- lm(y ~ group, data = dat)
  mse <- sum(residuals(fit)^2) / df.residual(fit)
  n_per_group <- table(dat$group)

  p_values <- numeric(n_pairs)
  for (i in 1:n_pairs) {
    g1 <- pairs[i, 1]
    g2 <- pairs[i, 2]

    diff <- abs(means[g1] - means[g2])
    se <- sqrt(mse * (1/n_per_group[g1] + 1/n_per_group[g2]))
    t_stat <- diff / se
    p_values[i] <- 2 * pt(-abs(t_stat), df = df.residual(fit))
  }

  # BH adjustment
  p_adjusted <- p.adjust(p_values, method = "BH")
  rejections <- p_adjusted < alpha

  list(
    any_reject = any(rejections),
    n_rejections = sum(rejections),
    p_values = p_adjusted,
    rejections = rejections
  )
}

## Scheffé method for all pairwise contrasts
scheffe_posthoc <- function(dat, alpha = 0.05) {
  k <- nlevels(dat$group)
  means <- tapply(dat$y, dat$group, mean)
  lev <- names(means)

  # All pairwise comparisons
  pairs <- t(combn(lev, 2))
  n_pairs <- nrow(pairs)

  # Pooled error from full ANOVA
  fit <- lm(y ~ group, data = dat)
  mse <- sum(residuals(fit)^2) / df.residual(fit)
  n_per_group <- table(dat$group)
  df_error <- df.residual(fit)

  # Scheffé critical value: sqrt((k-1) * F_{k-1, df_error, alpha})
  scheffe_crit <- sqrt((k - 1) * qf(1 - alpha, k - 1, df_error))

  rejections <- logical(n_pairs)
  for (i in 1:n_pairs) {
    g1 <- pairs[i, 1]
    g2 <- pairs[i, 2]

    diff <- abs(means[g1] - means[g2])
    se <- sqrt(mse * (1/n_per_group[g1] + 1/n_per_group[g2]))
    t_stat <- diff / se

    # Reject if |t| > Scheffé critical value
    rejections[i] <- t_stat > scheffe_crit
  }

  list(
    any_reject = any(rejections),
    n_rejections = sum(rejections),
    rejections = rejections
  )
}

## ------------------------------------------------------------
## 2. Wrapper for fusion cascade with calibrated c
## ------------------------------------------------------------

fusion_posthoc <- function(dat, c_factor, alpha = 0.05, return_partition = FALSE) {
  k <- nlevels(dat$group)
  m <- k - 1
  alpha_seq <- c_factor * alpha / (m:1)

  result <- cascade_fusions(dat, alpha_seq = alpha_seq)

  # Reconstruct the partition - always needed for counting differences
  # Start with each group in its own cluster
  partition <- as.list(1:k)

  if (!is.null(result$tests) && nrow(result$tests) > 0) {
    n_tests <- nrow(result$tests)
    # Apply fusions up to (but not including) the rejected one
    n_to_apply <- if (result$any_reject) n_tests - 1 else n_tests

    for (i in seq_len(n_to_apply)) {
      # g1 and g2 are strings like "1", "2", or "1_2" (compound labels)
      g1_str <- as.character(result$tests$g1[i])
      g2_str <- as.character(result$tests$g2[i])

      # Extract the original group numbers from compound labels
      # "1_2" -> c(1, 2)
      groups1 <- as.integer(strsplit(g1_str, "_")[[1]])
      groups2 <- as.integer(strsplit(g2_str, "_")[[1]])

      # Find which partition clusters contain any of these groups
      idx1 <- which(sapply(partition, function(cl) any(groups1 %in% cl)))
      idx2 <- which(sapply(partition, function(cl) any(groups2 %in% cl)))

      if (length(idx1) > 0 && length(idx2) > 0 && idx1[1] != idx2[1]) {
        # Merge the two clusters
        partition[[idx1[1]]] <- c(partition[[idx1[1]]], partition[[idx2[1]]])
        partition[[idx2[1]]] <- NULL
      }
    }
  }

  # Count pairwise differences from partition
  # Pairs in different clusters are "different"
  n_implied_differences <- 0
  n_clusters <- length(partition)
  if (n_clusters > 1) {
    for (i in 1:(n_clusters - 1)) {
      for (j in (i + 1):n_clusters) {
        # All pairs across clusters i and j are different
        n_implied_differences <- n_implied_differences +
          length(partition[[i]]) * length(partition[[j]])
      }
    }
  }

  out <- list(
    any_reject = result$any_reject,
    n_rejections = n_implied_differences,  # number of implied pairwise differences
    tests = result$tests
  )

  # If detailed partition info requested, add pairwise rejections vector
  if (return_partition) {
    # Convert partition to pairwise rejections
    pairs <- t(combn(k, 2))
    n_pairs <- nrow(pairs)
    rejections <- logical(n_pairs)

    for (i in 1:n_pairs) {
      g1 <- pairs[i, 1]
      g2 <- pairs[i, 2]
      # Check if g1 and g2 are in different clusters
      in_same_cluster <- any(sapply(partition, function(cl) g1 %in% cl && g2 %in% cl))
      rejections[i] <- !in_same_cluster
    }

    out$partition <- partition
    out$rejections <- rejections
  }

  out
}

## ------------------------------------------------------------
## 3. Power simulation scenarios
## ------------------------------------------------------------

## Scenario 1: Monotonic pattern (ordered means)
scenario_monotonic <- function(k = 5, effect_size = 0.5, sigma = 1) {
  # Means: 0, delta, 2*delta, ..., (k-1)*delta
  delta <- effect_size * sigma
  mu <- (0:(k-1)) * delta
  mu
}

## Scenario 2: Two-group pattern (k1 groups at 0, k2 groups at delta)
scenario_two_groups <- function(k = 6, k1 = 3, effect_size = 1.0, sigma = 1) {
  delta <- effect_size * sigma
  mu <- c(rep(0, k1), rep(delta, k - k1))
  mu
}

## Scenario 3: Three-group pattern
scenario_three_groups <- function(k = 9, k1 = 3, k2 = 3,
                                   effect_size1 = 0.5, effect_size2 = 1.0,
                                   sigma = 1) {
  delta1 <- effect_size1 * sigma
  delta2 <- effect_size2 * sigma
  mu <- c(rep(0, k1), rep(delta1, k2), rep(delta2, k - k1 - k2))
  mu
}

## Scenario 4: One outlier group
scenario_outlier <- function(k = 5, effect_size = 2.0, sigma = 1) {
  delta <- effect_size * sigma
  mu <- c(rep(0, k - 1), delta)
  mu
}

## ------------------------------------------------------------
## 4. Power assessment for one scenario
## ------------------------------------------------------------

assess_power <- function(mu, n = 20, sigma = 1, c_factor,
                         alpha = 0.05, R = 1000,
                         compute_accuracy = FALSE, tol = 0.01) {
  k <- length(mu)

  # Storage for results
  fusion_detect <- 0
  tukey_detect <- 0
  holm_detect <- 0
  bh_detect <- 0
  scheffe_detect <- 0

  fusion_n_rej <- numeric(R)
  tukey_n_rej <- numeric(R)
  holm_n_rej <- numeric(R)
  bh_n_rej <- numeric(R)
  scheffe_n_rej <- numeric(R)

  # If compute_accuracy, track true/false positives
  if (compute_accuracy) {
    # Determine true partition (groups with same mean, within tolerance)
    true_partition <- list()
    used <- rep(FALSE, k)

    for (i in 1:k) {
      if (!used[i]) {
        cluster <- which(abs(mu - mu[i]) < tol)
        true_partition[[length(true_partition) + 1]] <- cluster
        used[cluster] <- TRUE
      }
    }

    # Determine which pairs are truly different
    pairs <- t(combn(k, 2))
    n_pairs <- nrow(pairs)
    true_different <- logical(n_pairs)

    for (i in 1:n_pairs) {
      g1 <- pairs[i, 1]
      g2 <- pairs[i, 2]
      # Different if they're in different clusters of true partition
      in_same_cluster <- any(sapply(true_partition, function(cl) g1 %in% cl && g2 %in% cl))
      true_different[i] <- !in_same_cluster
    }

    n_true_diff <- sum(true_different)
    n_true_same <- sum(!true_different)

    # Storage for TP, FP, FN
    fusion_tp <- numeric(R)
    fusion_fp <- numeric(R)
    tukey_tp <- numeric(R)
    tukey_fp <- numeric(R)
    holm_tp <- numeric(R)
    holm_fp <- numeric(R)
    bh_tp <- numeric(R)
    bh_fp <- numeric(R)
    scheffe_tp <- numeric(R)
    scheffe_fp <- numeric(R)
  }

  for (r in 1:R) {
    dat <- simulate_data(mu, n = n, sigma = sigma)

    # Fusion cascade
    fusion_res <- fusion_posthoc(dat, c_factor = c_factor, alpha = alpha,
                                  return_partition = compute_accuracy)
    if (fusion_res$any_reject) fusion_detect <- fusion_detect + 1
    fusion_n_rej[r] <- fusion_res$n_rejections

    # Tukey HSD
    tukey_res <- tukey_posthoc(dat, alpha = alpha)
    if (tukey_res$any_reject) tukey_detect <- tukey_detect + 1
    tukey_n_rej[r] <- tukey_res$n_rejections

    # Holm
    holm_res <- holm_posthoc(dat, alpha = alpha)
    if (holm_res$any_reject) holm_detect <- holm_detect + 1
    holm_n_rej[r] <- holm_res$n_rejections

    # BH
    bh_res <- bh_posthoc(dat, alpha = alpha)
    if (bh_res$any_reject) bh_detect <- bh_detect + 1
    bh_n_rej[r] <- bh_res$n_rejections

    # Scheffé
    scheffe_res <- scheffe_posthoc(dat, alpha = alpha)
    if (scheffe_res$any_reject) scheffe_detect <- scheffe_detect + 1
    scheffe_n_rej[r] <- scheffe_res$n_rejections

    # Compute accuracy metrics if requested
    if (compute_accuracy) {
      # For fusion: use the reconstructed partition
      fusion_tp[r] <- sum(fusion_res$rejections & true_different)
      fusion_fp[r] <- sum(fusion_res$rejections & !true_different)

      # For pairwise methods: they directly return which pairs are rejected
      tukey_tp[r] <- sum(tukey_res$rejections & true_different)
      tukey_fp[r] <- sum(tukey_res$rejections & !true_different)

      holm_tp[r] <- sum(holm_res$rejections & true_different)
      holm_fp[r] <- sum(holm_res$rejections & !true_different)

      bh_tp[r] <- sum(bh_res$rejections & true_different)
      bh_fp[r] <- sum(bh_res$rejections & !true_different)

      scheffe_tp[r] <- sum(scheffe_res$rejections & true_different)
      scheffe_fp[r] <- sum(scheffe_res$rejections & !true_different)
    }
  }

  result <- data.frame(
    method = c("Fusion", "Tukey", "Holm", "BH", "Scheffe"),
    power = c(fusion_detect, tukey_detect, holm_detect, bh_detect, scheffe_detect) / R,
    mean_rejections = c(mean(fusion_n_rej), mean(tukey_n_rej),
                        mean(holm_n_rej), mean(bh_n_rej), mean(scheffe_n_rej))
  )

  if (compute_accuracy) {
    result$mean_tp <- c(mean(fusion_tp), mean(tukey_tp), mean(holm_tp), mean(bh_tp), mean(scheffe_tp))
    result$mean_fp <- c(mean(fusion_fp), mean(tukey_fp), mean(holm_fp), mean(bh_fp), mean(scheffe_fp))
    result$n_true_diff <- n_true_diff
    result$n_true_same <- n_true_same
  }

  result
}

## ------------------------------------------------------------
## 5. Accuracy of partition recovery (for fusion cascade)
## ------------------------------------------------------------

## Check if fusion cascade recovers the true partition
## This requires knowing the "true" partition structure
assess_partition_accuracy <- function(mu, n = 20, sigma = 1, c_factor,
                                      alpha = 0.05, R = 1000, tol = 0.01) {
  k <- length(mu)

  # Determine true partition (groups with same mean, within tolerance)
  true_partition <- list()
  used <- rep(FALSE, k)

  for (i in 1:k) {
    if (!used[i]) {
      cluster <- which(abs(mu - mu[i]) < tol)
      true_partition[[length(true_partition) + 1]] <- cluster
      used[cluster] <- TRUE
    }
  }

  n_true_clusters <- length(true_partition)

  # Run simulations
  correct_partition <- 0
  n_clusters_found <- numeric(R)

  for (r in 1:R) {
    dat <- simulate_data(mu, n = n, sigma = sigma)
    fusion_res <- fusion_posthoc(dat, c_factor = c_factor, alpha = alpha)

    # Count number of clusters in final model
    # If cascade stops at step j, we have k-j clusters
    if (fusion_res$any_reject && !is.null(fusion_res$tests)) {
      n_fusions <- nrow(fusion_res$tests) - 1  # last test was rejected
    } else {
      n_fusions <- k - 1  # all fused
    }

    n_clusters <- k - n_fusions
    n_clusters_found[r] <- n_clusters

    if (n_clusters == n_true_clusters) {
      correct_partition <- correct_partition + 1
    }
  }

  list(
    true_n_clusters = n_true_clusters,
    accuracy = correct_partition / R,
    mean_n_clusters = mean(n_clusters_found),
    sd_n_clusters = sd(n_clusters_found)
  )
}

## ------------------------------------------------------------
## 6. Example: Run power simulations for multiple scenarios
## ------------------------------------------------------------

run_power_study <- function(c_factor = 0.1, n = 20, sigma = 1,
                            alpha = 0.05, R = 1000) {

  cat("Running power simulations with n =", n, "R =", R, "\n\n")

  # Scenario 1: Monotonic
  cat("Scenario 1: Monotonic (k=5, effect=0.5)\n")
  mu1 <- scenario_monotonic(k = 5, effect_size = 0.5, sigma = sigma)
  power1 <- assess_power(mu1, n = n, sigma = sigma, c_factor = c_factor,
                         alpha = alpha, R = R)
  print(power1)
  cat("\n")

  # Scenario 2: Two groups
  cat("Scenario 2: Two groups (k=6, 3+3, effect=1.0)\n")
  mu2 <- scenario_two_groups(k = 6, k1 = 3, effect_size = 1.0, sigma = sigma)
  power2 <- assess_power(mu2, n = n, sigma = sigma, c_factor = c_factor,
                         alpha = alpha, R = R)
  print(power2)
  cat("\n")

  # Scenario 3: Three groups
  cat("Scenario 3: Three groups (k=9, 3+3+3, effects=0.5,1.0)\n")
  mu3 <- scenario_three_groups(k = 9, k1 = 3, k2 = 3,
                                effect_size1 = 0.5, effect_size2 = 1.0,
                                sigma = sigma)
  power3 <- assess_power(mu3, n = n, sigma = sigma, c_factor = c_factor,
                         alpha = alpha, R = R)
  print(power3)
  cat("\n")

  # Scenario 4: Outlier
  cat("Scenario 4: One outlier (k=5, effect=2.0)\n")
  mu4 <- scenario_outlier(k = 5, effect_size = 2.0, sigma = sigma)
  power4 <- assess_power(mu4, n = n, sigma = sigma, c_factor = c_factor,
                         alpha = alpha, R = R)
  print(power4)
  cat("\n")

  # Combine results
  results <- rbind(
    cbind(scenario = "Monotonic", power1),
    cbind(scenario = "TwoGroups", power2),
    cbind(scenario = "ThreeGroups", power3),
    cbind(scenario = "Outlier", power4)
  )

  results
}

## ------------------------------------------------------------
## 7. Example usage
## ------------------------------------------------------------

## For k=5, n=20, retrieve calibrated c from our data
## (You would load this from your CSV file)
# calib_data <- read.csv("nested_models_calibration_fusion.csv")
# c_val <- calib_data[calib_data$k == 5 & calib_data$n == 20, "c"]

## Example: use a typical value
set.seed(123)
# results <- run_power_study(c_factor = 0.1, n = 20, sigma = 1, R = 1000)
# write.csv(results, "power_comparison_results.csv", row.names = FALSE)

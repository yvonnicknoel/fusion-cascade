## ------------------------------------------------------------
## 0. Packages
## ------------------------------------------------------------
library(parallel)

## ------------------------------------------------------------
## 1. Data generation under H0 (all means equal)
## ------------------------------------------------------------
simulate_data <- function(mu, n = 20, sigma = 1) {
  k <- length(mu)
  group <- factor(rep(seq_len(k), each = n))
  y <- rnorm(k * n, mean = rep(mu, each = n), sd = sigma)
  data.frame(y = y, group = group)
}

## ------------------------------------------------------------
## 2. Fusion cascade (unconstrained version)
##    - alpha_seq: vector of thresholds (length k-1)
##      if NULL, a fixed alpha is used at each step
## ------------------------------------------------------------
cascade_fusions <- function(dat, alpha = 0.05, alpha_seq = NULL) {
  dat$group <- droplevels(dat$group)
  k <- nlevels(dat$group)

  # Full model (to fix the error term)
  M_full <- lm(y ~ group, data = dat)
  group_base <- dat$group

  # Threshold sequence
  if (is.null(alpha_seq)) {
    alpha_seq <- rep(alpha, k - 1)
  } else if (length(alpha_seq) < (k - 1)) {
    stop("alpha_seq too short.")
  }

  tests <- list()
  any_reject <- FALSE

  for (step in seq_len(k - 1)) {
    means <- tapply(dat$y, group_base, mean)
    lev   <- names(means)
    if (length(lev) <= 1L) break  # nothing left to fuse

    # Closest pair of means
    pairs <- t(combn(lev, 2))  # n_pairs x 2
    diff_means <- abs(means[pairs[, 1]] - means[pairs[, 2]])
    idx_min <- which.min(diff_means)
    g1 <- pairs[idx_min, 1]
    g2 <- pairs[idx_min, 2]

    # New fused factor
    new_group <- group_base
    levels(new_group)[levels(new_group) %in% c(g1, g2)] <- paste(g1, g2, sep = "_")
    new_group <- droplevels(new_group)

    M_base <- lm(y ~ group_base, data = dat)
    if (nlevels(new_group) > 1L) {
      M_cand <- lm(y ~ new_group, data = dat)
    } else {
      M_cand <- lm(y ~ 1, data = dat)  # null model
    }

    # 3-model ANOVA: M_cand, M_base, M_full
    tab    <- anova(M_cand, M_base, M_full)
    F_stat <- tab$"F"[2]
    p_val  <- tab$"Pr(>F)"[2]

    tests[[step]] <- data.frame(
      step      = step,
      g1        = g1,
      g2        = g2,
      diff_mean = diff_means[idx_min],
      F         = F_stat,
      p         = p_val,
      alpha     = alpha_seq[step]
    )

    if (!is.na(p_val) && p_val < alpha_seq[step]) {
      any_reject <- TRUE
      break
    } else {
      group_base <- new_group
    }
  }

  list(
    any_reject = any_reject,
    tests      = if (length(tests)) do.call(rbind, tests) else NULL
  )
}

## ------------------------------------------------------------
## 3. Calibration of the cascade for a given (k, n)
##    -> finds c such that FWER ~ target_alpha under H0
## ------------------------------------------------------------
calibrate_cascade <- function(k = 5, n = 20, sigma = 1,
                              target_alpha = 0.05,
                              R = 2000,
                              tol = 0.005,
                              max_iter = 25,
                              c_min = 1e-8,  # lower bound for c
                              c_max = 1      # upper bound for c
                              ) {
  mu_null <- rep(0, k)
  m <- k - 1

  # Helper function: FWER for a given c factor
  fwer_for_c <- function(c_factor) {
    alpha <- target_alpha
    alpha_seq <- c_factor * alpha / (m:1)
    count_reject <- 0L

    for (r in 1:R) {
      dat <- simulate_data(mu_null, n = n, sigma = sigma)
      cas <- cascade_fusions(dat, alpha_seq = alpha_seq)
      if (cas$any_reject) count_reject <- count_reject + 1L
    }

    count_reject / R
  }

  ## Binary search on log10(c) between log10(c_min) and log10(c_max)
  log_lower <- log10(c_min)
  log_upper <- log10(c_max)

  best_c    <- NA_real_
  best_fwer <- NA_real_

  for (iter in 1:max_iter) {
    log_mid  <- (log_lower + log_upper) / 2
    c_mid    <- 10^log_mid
    fwer_mid <- fwer_for_c(c_mid)
    # cat("iter", iter, "c =", c_mid, "FWER ~", fwer_mid, "\n")

    best_c    <- c_mid
    best_fwer <- fwer_mid

    if (abs(fwer_mid - target_alpha) < tol) {
      break
    }

    if (fwer_mid > target_alpha) {
      # Too liberal -> reduce c -> go down in log10(c)
      log_upper <- log_mid
    } else {
      # Too conservative -> increase c -> go up in log10(c)
      log_lower <- log_mid
    }
  }

  list(
    c_factor = best_c,
    fwer     = best_fwer
  )
}

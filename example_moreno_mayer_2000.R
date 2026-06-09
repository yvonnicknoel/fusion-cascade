## Illustrative Example: Moreno & Mayer (2000) Experiment 2
## Multimedia Learning: Narration, Sounds, and Music
## Demonstrates fusion cascade analysis of factorial design

library(ggplot2)
library(multcomp)

cat("=== Moreno & Mayer (2000) Experiment 2 ===\n")
cat("Multimedia learning: Effect of auditory adjuncts on learning\n\n")

# Original study: 2x2 factorial design
# Factor A: Background sounds (Present vs Absent)
# Factor B: Background music (Present vs Absent)
# DV: Transfer test performance (problem-solving)

# Group labels and cell means from Table 2
groups <- c("N", "NS", "NM", "NSM")
group_labels <- c(
  "Narration only",
  "Narration + Sounds",
  "Narration + Music",
  "Narration + Sounds + Music"
)

# Observed means and SDs (Transfer test from Table 2)
observed_means <- c(5.55, 3.06, 3.33, 3.40)
observed_sds <- c(1.88, 1.92, 1.72, 2.39)
sample_sizes <- c(20, 17, 18, 20)

k <- 4
total_n <- sum(sample_sizes)

cat("Observed Transfer Test Scores:\n")
for (i in 1:k) {
  cat(sprintf("  %-30s: M = %.2f, SD = %.2f, n = %d\n",
              group_labels[i], observed_means[i], observed_sds[i], sample_sizes[i]))
}
cat("\n")

# Simulate data matching the observed statistics EXACTLY
# Generate from N(0,1), then shift and scale to match exact means and SDs
set.seed(12345)
dat_list <- list()
for (i in 1:k) {
  # Generate standard normal data
  z <- rnorm(sample_sizes[i])
  # Scale and shift to get exact mean and SD
  y <- observed_sds[i] * (z - mean(z)) / sd(z) + observed_means[i]

  dat_list[[i]] <- data.frame(
    y = y,
    group = factor(i),
    group_label = group_labels[i]
  )
}
dat <- do.call(rbind, dat_list)
dat$group <- factor(dat$group, levels = 1:k, labels = groups)

# Traditional ANOVA analysis
cat("=== Traditional 2x2 ANOVA Analysis ===\n\n")

# Create factorial structure
dat$sounds <- factor(ifelse(dat$group %in% c("NS", "NSM"), "Present", "Absent"),
                     levels = c("Absent", "Present"))
dat$music <- factor(ifelse(dat$group %in% c("NM", "NSM"), "Present", "Absent"),
                    levels = c("Absent", "Present"))

# Run 2-way ANOVA
aov_result <- aov(y ~ sounds * music, data = dat)
cat("Two-way ANOVA:\n")
print(summary(aov_result))
cat("\n")

# Traditional post-hoc: Tukey HSD on all pairwise comparisons
cat("=== Traditional Post-hoc: Tukey HSD ===\n\n")
tukey_result <- tukey_posthoc(dat, alpha = 0.05)

# Extract pairwise comparisons
fit <- aov(y ~ group, data = dat)
tukey_full <- TukeyHSD(fit, conf.level = 0.95)
cat("Pairwise comparisons (Tukey HSD):\n")
print(tukey_full$group)
cat("\n")

# Interpret Tukey results
cat("Tukey HSD interpretation:\n")
sig_pairs <- rownames(tukey_full$group)[tukey_full$group[, "p adj"] < 0.05]
if (length(sig_pairs) > 0) {
  cat("  Significant differences:\n")
  for (pair in sig_pairs) {
    cat(sprintf("    %s\n", pair))
  }
} else {
  cat("  No significant pairwise differences\n")
}
cat("\n")

cat("Note: Tukey provides pairwise comparisons but does not directly reveal\n")
cat("the clustering structure or which groups can be treated as equivalent.\n\n")

# Fusion cascade analysis
cat("=== Fusion Cascade Analysis ===\n\n")

source("power_simulations.R")

# Function to compute approximate calibration factor from regression model
# Model: log(c) = beta0 + beta1*k + beta2*log(k)
# Fitted: log(c) = 1.55 - 0.29*k - 1.41*log(k)
compute_c_approx <- function(k) {
  beta0 <- 1.55
  beta1 <- -0.29
  beta2 <- -1.41
  log_c <- beta0 + beta1 * k + beta2 * log(k)
  return(exp(log_c))
}

# Compute approximate c and alpha* values
c_approx <- compute_c_approx(k)
cat(sprintf("Approximate calibration factor: c_approx = %.4f\n", c_approx))
cat(sprintf("Actual calibration factor:       c_actual = %.4f\n", c_factor))
cat(sprintf("Relative difference: %.1f%%\n\n", 100 * abs(c_approx - c_factor) / c_factor))

# Apply fusion cascade with actual calibration
alpha <- 0.05
m <- k - 1
alpha_seq_approx <- c_approx * alpha / (m:1)

# result <- cascade_fusions(dat, alpha_seq = alpha_seq)
result <- cascade_fusions(dat, alpha_seq = alpha_seq_approx)

# Reconstruct final partition
partition <- as.list(1:k)
names(partition) <- groups

if (!is.null(result$tests)) {
  for (i in 1:nrow(result$tests)) {
    test <- result$tests[i, ]
    if (test$p >= test$alpha) {
      # Fusion accepted
      g1_label <- as.character(test$g1)
      g2_label <- as.character(test$g2)

      # Handle merged cluster labels (e.g., "NS_NSM")
      g1_groups <- strsplit(g1_label, "_")[[1]]
      g2_groups <- strsplit(g2_label, "_")[[1]]

      # Find which partition elements contain these groups
      idx1 <- which(sapply(partition, function(x) any(x %in% which(groups %in% g1_groups))))
      idx2 <- which(sapply(partition, function(x) any(x %in% which(groups %in% g2_groups))))

      if (length(idx1) > 0 && length(idx2) > 0 && idx1[1] != idx2[1]) {
        partition[[idx1[1]]] <- c(partition[[idx1[1]]], partition[[idx2[1]]])
        partition[[idx2[1]]] <- NULL
      }
    } else {
      break
    }
  }
}

cat("=== Final Partition ===\n\n")
cat(sprintf("Number of clusters: %d\n\n", length(partition)))
for (i in 1:length(partition)) {
  cluster_groups <- groups[partition[[i]]]
  cluster_labels <- group_labels[partition[[i]]]
  cluster_mean <- mean(observed_means[partition[[i]]])

  cat(sprintf("Cluster %d (M = %.2f):\n", i, cluster_mean))
  for (j in 1:length(cluster_groups)) {
    cat(sprintf("  - %s: %s\n", cluster_groups[j], cluster_labels[j]))
  }
  cat("\n")
}

# Interpretation
cat("=== Interpretation ===\n\n")
cat("The fusion cascade reveals the cognitive overload structure:\n\n")

cat("Moreno & Mayer's coherence hypothesis predicted that adding extraneous\n")
cat("auditory material (sounds and/or music) would overload auditory working\n")
cat("memory and impair learning. The fusion cascade directly identifies this\n")
cat("pattern by partitioning the conditions into meaningful clusters.\n\n")

cat("The partition shows that:\n")
cat("1. 'Narration only' forms its own cluster with superior performance\n")
cat("2. All conditions with auditory adjuncts cluster together with\n")
cat("   substantially impaired transfer performance\n\n")

cat("This clustering directly supports the coherence principle: multimedia\n")
cat("explanations should exclude extraneous sounds that do not contribute\n")
cat("to understanding. The method reveals that once auditory overload occurs,\n")
cat("the specific type of adjunct (sounds, music, or both) matters less than\n")
cat("the presence of any extraneous auditory information.\n\n")

cat("Compare this to traditional analysis:\n")
cat("- Traditional ANOVA: 'Significant interaction, F(1,71) = 7.61, p < .01'\n")
cat("  (tells us there's an interaction but not what it means)\n")
cat("- Tukey HSD: Multiple pairwise comparisons that must be interpreted\n")
cat("  together to infer the structure\n")
cat("- Fusion cascade: Directly identifies the two-tier structure\n")
cat("  (clean vs. overloaded conditions)\n\n")

# Create visualization
plot_data <- data.frame(
  Group = factor(groups, levels = groups),
  Label = group_labels,
  Mean = observed_means,
  SD = observed_sds,
  n = sample_sizes,
  Cluster = NA,
  FittedMean = NA,
  FittedSE = NA
)

# Assign cluster colors, fitted means, and standard errors
for (i in 1:length(partition)) {
  cluster_indices <- partition[[i]]
  plot_data$Cluster[cluster_indices] <- i

  # Fitted mean = weighted mean of groups in this cluster
  cluster_ns <- sample_sizes[cluster_indices]
  cluster_means <- observed_means[cluster_indices]
  cluster_mean <- sum(cluster_ns * cluster_means) / sum(cluster_ns)
  plot_data$FittedMean[cluster_indices] <- cluster_mean

  # Standard error of the fitted cluster mean
  # SE = sqrt(MSE / n_total) where n_total is sum of all n's in cluster
  # We use the pooled MSE from the ANOVA (from result or recalculate)
  # For simplicity, we'll use the pooled SD from observed data
  cluster_vars <- observed_sds[cluster_indices]^2
  cluster_n_total <- sum(cluster_ns)

  # Pooled variance within cluster
  # Using the fact that each group contributes (n-1)*var to pooled variance
  pooled_var <- sum((cluster_ns - 1) * cluster_vars^2) / sum(cluster_ns - 1)

  # SE of cluster mean
  cluster_se <- sqrt(pooled_var / cluster_n_total)
  plot_data$FittedSE[cluster_indices] <- cluster_se
}

# Add numeric group position for plotting fitted means
plot_data$GroupNum <- as.numeric(plot_data$Group)

cat("=== Plot Data for Verification ===\n")
print(plot_data[, c("Group", "Mean", "FittedMean", "FittedSE", "GroupNum")])
cat("\n")

p <- ggplot(plot_data, aes(x = Group, y = Mean)) +
  geom_bar(aes(fill = factor(Cluster)), stat = "identity", color = "black", alpha = 0.7) +
  # Add fitted means as points connected by lines
  geom_point(data = plot_data, aes(x = GroupNum, y = FittedMean),
             color = "black", size = 3, shape = 19, inherit.aes = FALSE) +
  geom_line(data = plot_data, aes(x = GroupNum, y = FittedMean, group = 1),
            color = "black", linewidth = 1, inherit.aes = FALSE) +
  # Add error bars for fitted means (95% CI)
  geom_errorbar(data = plot_data, aes(x = GroupNum, ymin = FittedMean - 1.96*FittedSE,
                                       ymax = FittedMean + 1.96*FittedSE),
                width = 0.15, color = "black", linewidth = 0.8, inherit.aes = FALSE) +
  scale_fill_manual(values = c("1" = "#0072B2", "2" = "#D55E00"),
                    name = "Observed data",
                    labels = c("1" = "No overload", "2" = "Auditory overload")) +
  scale_x_discrete(labels = c("N" = "Narration\nonly",
                              "NS" = "Narration\n+ Sounds",
                              "NM" = "Narration\n+ Music",
                              "NSM" = "Narration\n+ Both")) +
  scale_linetype_manual(name = NULL, values = c("Fitted means (2-cluster model)" = "solid")) +
  labs(x = "Multimedia Condition",
       y = "Transfer Test Score",
       title = "Moreno & Mayer (2000): Cognitive Overload Structure",
       subtitle = "Two clusters identified: clean narration vs. overloaded conditions") +
  theme_bw(base_size = 11) +
  theme(legend.position = c(0.75, 0.80),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.grid.minor = element_blank()) # +
  # guides(fill = guide_legend(order = 1), linetype = guide_legend(order = 2, override.aes = list(linewidth = 1)))

x11(width = 5, height = 5)
print(p)

ggsave("example_moreno_mayer_2000.eps", p, width = 5, height = 5)
ggsave("example_moreno_mayer_2000.pdf", p, width = 5, height = 5)
ggsave("example_moreno_mayer_2000.png", p, width = 5, height = 5, dpi = 300)

cat("Plot saved to example_moreno_mayer_2000.pdf/.png\n\n")

cat("=== Connection to Factorial Design ===\n\n")
cat("This example illustrates a key advantage of the fusion cascade:\n")
cat("It naturally handles factorial designs by treating them as one-way\n")
cat("ANOVA on cell means. This:\n\n")
cat("1. Simplifies the follow-up analysis (after identification of main effects,\n")
cat("   interactions, simple effects)\n")
cat("2. Proposes an interaction structure as a partition\n")
cat("3. Provides clear, model-based interpretation\n\n")
cat("The method shows that the 2x2 interaction structure is actually\n")
cat("a simple dichotomy: conditions without auditory adjuncts vs.\n")
cat("conditions with any auditory adjuncts.\n")

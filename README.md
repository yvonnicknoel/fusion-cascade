# Fusion Cascade

**A sequential fusion procedure for post-hoc comparisons in fixed-effect ANOVA, with Monte-Carlo-calibrated error control.**

This repository accompanies the article:

> Noel, Y. (2026). *Post-hoc comparisons in fixed-effect ANOVA: A sequential fusion approach with calibrated error control.* **Psychological Methods.**

It contains the reference R implementation of the procedure, the Monte Carlo calibration data, and all scripts needed to reproduce the simulations and the illustrative example reported in the paper.

---

## What the method does

Standard post-hoc procedures (Tukey HSD, Holm–Bonferroni, Benjamini–Hochberg, Scheffé) treat the analysis as a collection of independent pairwise decisions. This can produce **non-transitive** verdicts — group A is judged equal to B, B equal to C, yet A differs from C — which correspond to no coherent grouping of the conditions.

The **fusion cascade** reframes post-hoc analysis as a model-selection problem. Starting from the saturated ANOVA model (one mean per group), it:

1. **Fuses** the pair of groups (or clusters) with the closest means;
2. **Tests** each fusion against the saturated-model error term via a likelihood-ratio / *F*-test;
3. **Stops** at the first fusion that significantly worsens fit, using a calibrated sequence of critical values.

The output is always a **single, transitive partition** of the groups — directly interpretable and representable as a dendrogram — while keeping the family-wise error rate (FWER) at the nominal level.

The calibrated critical values follow

```
αⱼ*(k, n) = c(k, n) · α / (m − j + 1),   j = 1, …, m = k − 1
```

with a simple closed-form approximation for the calibration factor:

```
c(k) ≈ 4.7 · exp(−0.29 k) · k^(−1.41)
```

> The cascade is an **exploratory** tool, offered as an improvement on other sequential post-hoc procedures. It complements — rather than replaces — planned contrasts when theory yields specific predictions.

---

## Repository contents

### Core implementation
| File | Description |
|------|-------------|
| `nested_models.R` | Core functions: `simulate_data()`, `cascade_fusions()` (the sequential procedure), and `calibrate_cascade()` (Monte Carlo binary search for *c(k, n)*). |
| `power_simulations.R` | Comparison framework implementing Tukey HSD, Holm–Bonferroni, Benjamini–Hochberg, and Scheffé alongside the cascade. Sources `nested_models.R`. |

### Simulation scripts (reproduce the paper's tables)
| File | Produces |
|------|----------|
| `run_fwer_null.R` | FWER under the global null for all methods (Table 3). |
| `run_power_comparisons-capprox.R` | Full power comparison across scenarios and sample sizes (Tables 4–5). |
| `assess_accuracy.R` | True- and false-positive rates across scenarios (Table 6). |
| `run_interpretational_coherence_capprox.R` | Partition recovery and non-transitivity rates vs. Tukey HSD (Table 7). |

### Illustrative example
| File | Description |
|------|-------------|
| `example_moreno_mayer_2000.R` | Reanalysis of Moreno & Mayer (2000, Exp. 2), a 2×2 multimedia-learning design (Section 6). |

### Data files
| File | Description |
|------|-------------|
| `nested_models_calibration_fusion.csv` | Calibrated *c(k, n)* values and empirical FWER for *k* = 3–50, *n* = 5–50 (Table 2). |
| `power_comparison_results_capprox.csv` | Power results across all scenarios and methods (Tables 4–5). |
| `fwer_null_results.csv` | Empirical FWER under the global null (Table 3). |

---

## Requirements

R (≥ 4.0) with the packages:

```r
install.packages(c("parallel", "multcomp", "ggplot2"))
```

---

## Quick start

Apply the fusion cascade to your own balanced one-way design:

```r
source("nested_models.R")

# Your data: a data.frame with a numeric column 'y' and a factor 'group'
dat <- data.frame(y = ..., group = factor(...))

# Calibration factor for your number of groups
k        <- nlevels(dat$group)
c_approx <- 4.7 * exp(-0.29 * k) * k^(-1.41)

# Build the decreasing threshold sequence (target FWER = 0.05)
m         <- k - 1
alpha_seq <- c_approx * 0.05 / (m:1)

# Run the cascade
result <- cascade_fusions(dat, alpha_seq = alpha_seq)
print(result$tests)      # step-by-step fusions, F, p, and threshold
result$any_reject        # TRUE if the cascade stopped at a significant fusion
```

For an exact (rather than approximated) threshold for a specific *(k, n)*, calibrate directly:

```r
cal <- calibrate_cascade(k = 5, n = 20, target_alpha = 0.05, R = 2000)
cal$c_factor   # calibrated c
cal$fwer       # achieved empirical FWER
```

---

## Reproducing the paper's results

Each `run_*.R` script is self-contained and sources its own dependencies:

```r
source("run_power_comparisons-capprox.R")
```

> ⚠️ The simulation scripts use parallel processing and may take several hours depending on hardware.

---

## Citing

If you use this software, please cite both the paper and the software:

```bibtex
@article{Noel2026FusionCascade,
  author  = {Noel, Yvonnick},
  title   = {Post-hoc Comparisons in Fixed-Effect {ANOVA}:
             A Sequential Fusion Approach with Calibrated Error Control},
  journal = {Psychological Methods},
  year    = {2026}
}

@software{NoelFusionCascade2026,
  author    = {Noel, Yvonnick},
  title     = {Fusion Cascade: An {R} Implementation of the Sequential
               Fusion Procedure for Post-Hoc {ANOVA} Comparisons},
  year      = {2026},
  version   = {1.0},
  publisher = {GitHub},
  url       = {https://github.com/yvonnicknoel/fusion-cascade}
}
```

---

## License

Released under the [MIT License](LICENSE). You may reuse, modify, and redistribute
the code with attribution.

## Contact

Yvonnick Noel — Université Rennes 2, France
ORCID: [0000-0002-3363-1365](https://orcid.org/0000-0002-3363-1365)

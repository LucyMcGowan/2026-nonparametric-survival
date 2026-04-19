

# Nonparametric Alternative for Summarizing Survival Treatment Effects

**D’Agostino McGowan, Rigdon, Li, and Small**

This repository contains code and data to reproduce the analyses in the
paper. The two main components are the simulation study and the real
data example using the rhDNase trial.

## Reproducing the Simulations

The simulation study compares the proposed randomization-based
confidence intervals for the multiplicative factor against Weibull AFT
model-based intervals across two data-generating processes (Weibull and
log-logistic), two sample sizes (50, 200), three censoring rates (20%,
50%, 80%), and two effect sizes (1.0 representing the null and 0.8).

To reproduce the simulation results and figures:

``` r
# Run the full simulation (1,000 replicates per scenario)
source("sim/run_simulations.R")

# Reproduce Figures 1–3
source("sim/plot_simulations.R")
```

- `sim/run_simulations.R` — generates all simulation results and saves
  them to `sim/results/`
- `sim/plot_simulations.R` — reads the saved results and produces
  `fig_type1.png`, `fig_coverage.png`, and `fig_width.png`

## Real Data Example

We illustrate the proposed confidence intervals using the rhDNase trial,
a randomized double-blind study of recombinant human deoxyribonuclease I
(rhDNase) for cystic fibrosis. The primary endpoint is time to first
pulmonary exacerbation. The dataset is available in the `survival` R
package.

### Setup

``` r
library(survival)

dat_first <- subset(rhDNase, !duplicated(id))
dat <- data.frame(
  Y = ifelse(!is.na(dat_first$ivstart), dat_first$ivstart,
             as.numeric(dat_first$end.dt - dat_first$entry.dt)),
  Z = dat_first$trt,
  event = as.integer(!is.na(dat_first$ivstart))
)
dat <- dat[dat$Y > 0, ]  # exclude 6 subjects infected at enrollment
```

Of the 641 remaining subjects, 323 received placebo and 318 received
rhDNase, with 138 and 103 events respectively. The log-rank test gives p
= 0.0056.

### Additive Shift CI

The function below inverts the log-rank test over a grid of candidate
additive shift values $c$ to produce a confidence interval for the
number of days by which treatment delays the event.

``` r
ci_additive <- function(data, lower, upper, alpha = 0.05, ngrid = 1001) {
  grid_c <- seq(lower, upper, length.out = ngrid)
  pval <- vapply(grid_c, function(c) {
    d <- data
    d$Y[d$Z == 0] <- d$Y[d$Z == 0] + c
    survdiff(Surv(Y, event) ~ Z, data = d, rho = 0)$pvalue
  }, numeric(1))

  data.frame(
    est   = grid_c[which.max(pval)],
    lower = min(grid_c[pval > alpha]),
    upper = max(grid_c[pval > alpha])
  )
}

ci_additive(dat, -80, 80)
```

        est lower upper
    1 50.08 17.12 72.96

The 95% CI is **\[17.1, 73.0\] days**, with a point estimate of ~50
days. rhDNase is estimated to have delayed time to first exacerbation by
roughly two to ten weeks.

### Multiplicative Factor CI

The function below inverts the log-rank test over a grid of candidate
multiplicative factors (searched on the log scale) to produce a
confidence interval for the time-acceleration factor.

``` r
ci_multiplicative <- function(data, log_lower, log_upper, alpha = 0.05, ngrid = 1001) {
  grid_rho <- exp(seq(log_lower, log_upper, length.out = ngrid))
  pval <- vapply(grid_rho, function(rho) {
    d <- data
    d$Y[d$Z == 0] <- d$Y[d$Z == 0] / rho
    survdiff(Surv(Y, event) ~ Z, data = d, rho = 0)$pvalue
  }, numeric(1))

  data.frame(
    est   = grid_rho[which.max(pval)],
    lower = min(grid_rho[pval > alpha]),
    upper = max(grid_rho[pval > alpha])
  )
}

ci_multiplicative(dat, -3, 3)
```

            est     lower    upper
    1 0.6453258 0.5294058 0.876341

The point estimate is 0.645 with 95% CI **\[0.529, 0.876\]**.

### Kaplan–Meier Plot with Transformed Control Curves

To assess whether a constant-effect summary is adequate, we overlay two
transformations of the control arm’s Kaplan–Meier curve: one shifted
right by $\hat{c} = 50$ days (additive) and one with its time axis
stretched by $1/\hat\rho \approx 1.55$ (multiplicative).

``` r
library(ggsurvfit)
```

    Loading required package: ggplot2

``` r
library(ggplot2)

km_fit     <- survfit2(Surv(Y, event) ~ Z, data = dat)
km_control <- survfit2(Surv(Y, event) ~ 1, data = subset(dat, Z == 0)) |>
  tidy_survfit()

km_mult <- km_control |> transform(time = time / 0.645)
km_add  <- km_control |> transform(time = time + 50)

km_fit |>
  ggsurvfit() +
  geom_step(
    data = km_mult, aes(x = time, y = estimate),
    linetype = "dashed", inherit.aes = FALSE
  ) +
  geom_step(
    data = km_add, aes(x = time, y = estimate),
    linetype = "dotted", inherit.aes = FALSE
  ) +
  scale_ggsurvfit(x_scales = list(limits = c(0, 189)))
```

    Warning: Removed 32 rows containing missing values or values outside the scale range
    (`geom_step()`).

    Warning: Removed 26 rows containing missing values or values outside the scale range
    (`geom_step()`).

![Kaplan–Meier survival curves for the rhDNase trial. Solid lines show
the control (blue) and rhDNase (orange) arms. The dashed line is the
control curve with its time axis stretched by $1/\hat\rho \approx 1.55$
(multiplicative); the dotted line is the control curve shifted right by
$\hat{c} = 50$ days (additive). The multiplicative transformation tracks
the treatment arm more closely, suggesting a time-acceleration factor is
the more adequate
summary.](Readme_files/figure-commonmark/unnamed-chunk-4-1.png)

## Shiny App

An interactive Shiny application for computing both confidence intervals
is available at: <https://lucy.shinyapps.io/survival_ci>

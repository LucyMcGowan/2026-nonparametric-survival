library(tidyverse)

load("sim/results/results.rda")

long <- results |>
  pivot_longer(
    cols = c(np_coverage, aft_coverage, np_width, aft_width),
    names_to  = c("method", "metric"),
    names_sep = "_",
    values_to = "value"
  ) |>
  mutate(
    method = recode(method, np = "Nonparametric\nRandomization-based", aft = "Weibull AFT"),
    dgp    = recode(dgp, weibull = "Weibull", loglogistic = "Log-logistic")
  )

# Plot 1: Type I error (rho = 1, plot 1 - coverage)
p1 <- long |>
  filter(metric == "coverage", rho == 1) |>
  mutate(type1 = 1 - value) |>
  ggplot(aes(
    x = target_cens,
    y = type1,
    color = method,
    group = method
  )) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dotted",
    color = "grey50",
    linewidth = 0.6
  ) +
  scale_color_manual(values = c("cornflower blue", "orange")) + 
  scale_x_continuous(breaks = c(0.2, 0.5, 0.8)) + 
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2,
             fill = "white",
             stroke = 1.2) +
  facet_grid(dgp ~ n, labeller = label_both) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Censoring rate", y = "Type I error", color = NULL) +
  theme_minimal() + 
  theme(legend.position = "bottom")


# Plot 2: Coverage (rho = 0.8, x = cens_rate, facet by n)
p2 <- long |>
  filter(metric == "coverage", rho == 0.8) |>
  ggplot(aes(
    x = target_cens,
    y = value,
    color = method,
    group = method
  )) +
  geom_hline(
    yintercept = 0.95,
    linetype = "dotted",
    color = "grey50",
    linewidth = 0.6
  ) +
  scale_color_manual(values = c("cornflower blue", "orange")) + 
  scale_x_continuous(breaks = c(0.2, 0.5, 0.8)) + 
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2,
             fill = "white",
             stroke = 1.2) +
  facet_grid(dgp ~ n, labeller = label_both) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Censoring rate", y = "Coverage", color = NULL) +
  theme_minimal() + 
  theme(legend.position = "bottom")


# Plot 3: CI width (rho = 0.8, x = cens_rate, facet by n)
p3 <- long |>
  filter(metric == "width", rho == 0.8) |>
  ggplot(aes(
    x = target_cens,
    y = value,
    color = method,
    group = method
  )) +
  scale_color_manual(values = c("cornflower blue", "orange")) + 
  scale_x_continuous(breaks = c(0.2, 0.5, 0.8)) + 
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2,
             fill = "white",
             stroke = 1.2) +
  facet_grid(dgp ~ n, labeller = label_both) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Censoring rate", y = "CI width", color = NULL) +
  theme_minimal() + 
  theme(legend.position = "bottom")

p1
ggsave(file = "fig_type1.png")
p2
ggsave(file = "fig_coverage.png")
p3
ggsave(file = "fig_width.png")

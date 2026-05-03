library(tidyverse)
load("sim/results/results.rda")

long <- results |>
  pivot_longer(
    cols = c(np_coverage, aft_coverage, np_width, aft_width,),
    names_to  = c("method", "metric"),
    names_sep = "_",
    values_to = "value"
  ) |>
  mutate(
    method = recode(method, np = "Nonparametric\nRandomization-based", aft = "Weibull AFT"),
    dgp    = recode(dgp, weibull = "Weibull", loglogistic = "Log-logistic")
  )

plot_cols <- c("#0072B2", "#E69F00")

pub_theme <- theme_bw(base_size = 13) +
  theme(
    legend.position       = "bottom",
    legend.key.width      = unit(1.8, "cm"),
    panel.grid.minor      = element_blank(),
    panel.grid.major.x    = element_blank(),
    panel.grid.major.y    = element_line(color = "grey88", linewidth = 0.4),
    strip.background      = element_rect(fill = "white", color = "white"),
    strip.text            = element_text(face = "bold")
  )

pub_scales <- list(
  scale_color_manual(values = plot_cols),
  scale_linetype_manual(values = c("solid", "22")),
  scale_shape_manual(values = c(16, 21)),
  scale_x_continuous(breaks = c(0.2, 0.5, 0.8),
                     labels = scales::percent_format(accuracy = 1))
)

facet_labels <- labeller(
  n   = function(x) paste0("n = ", x),
  dgp = identity
)

# Plot 1: Type I error
p1 <- long |>
  filter(metric == "coverage", rho == 1) |>
  mutate(type1 = 1 - value) |>
  ggplot(aes(x = target_cens, y = type1,
             color = method, linetype = method,
             shape = method, group = method)) +
  geom_hline(yintercept = 0.05, linetype = "dotted") +
  pub_scales +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.4, fill = "white", stroke = 1.2) +
  facet_grid(dgp ~ n, labeller = facet_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = c(0.05, 0.1, 0.15)) +
  labs(x = "Censoring rate", y = "Type I error",
       color = NULL, linetype = NULL, shape = NULL) +
  pub_theme

# Plot 2: Coverage
p2 <- long |>
  filter(metric == "coverage", rho == 1.25) |>
  ggplot(aes(x = target_cens, y = value,
             color = method, linetype = method,
             shape = method, group = method)) +
  geom_hline(yintercept = 0.95, linetype = "dotted") +
  pub_scales +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.4, fill = "white", stroke = 1.2) +
  facet_grid(dgp ~ n, labeller = facet_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = c(0.86, 0.89, 0.92, 0.95, 0.98)) +
  labs(x = "Censoring rate", y = "Coverage probability",
       color = NULL, linetype = NULL, shape = NULL) +
  pub_theme

# Plot 3: CI width
p3 <- long |>
  filter(metric == "width", rho == 1.25) |>
  ggplot(aes(x = target_cens, y = value,
             color = method, linetype = method,
             shape = method, group = method)) +
  pub_scales +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.4, fill = "white", stroke = 1.2) +
  facet_grid(dgp ~ n, labeller = facet_labels) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Censoring rate", y = "Median CI width",
       color = NULL, linetype = NULL, shape = NULL) +
  pub_theme

ggsave("fig_type1.png",    p1, width = 6.5, height = 4.5, dpi = 300)
ggsave("fig_coverage.png", p2, width = 6.5, height = 4.5, dpi = 300)
ggsave("fig_width.png",    p3, width = 6.5, height = 4.5, dpi = 300)
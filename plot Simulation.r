
library(dplyr)
library(ggplot2)
library(patchwork)

setwd("Your Directory")

prev <- "0.10"
h2   <- "0.50"

AM <- read.table(paste0("simulation_AM1000_prev_", prev, "_h2_", h2, ".txt"), header = TRUE)
RM <- read.table(paste0("simulation_RM1000_prev_", prev, "_h2_", h2, ".txt"), header = TRUE)

RM$group <- "Random Mating CIF"
AM$group <- "Assortative Mating CIF"
data <- rbind(AM, RM)

# -----------------------------
# Colors
# -----------------------------
group_cols <- c(
  "Assortative Mating CIF" = "#556B2F",
  "Random Mating CIF"      = "#FFB90F"
)

# -----------------------------
# P-value label (mean test)
# -----------------------------
fmt_p <- function(p) {
  if (is.na(p)) return("t-test p = NA")
  if (p < 1e-4) return("t-test p < 1e-4")
  paste0("t-test p = ", signif(p, 3))
}

# -----------------------------
# Means (for vlines)
# -----------------------------
means_prev <- data %>%
  group_by(group) %>%
  summarise(mu = mean(prev, na.rm = TRUE), .groups = "drop")

means_rho <- data %>%
  group_by(group) %>%
  summarise(mu = mean(rho, na.rm = TRUE), .groups = "drop")

# -----------------------------
# t-tests
# -----------------------------
p_prev <- t.test(prev ~ group, data = data)$p.value
p_rho  <- t.test(rho  ~ group, data = data)$p.value

# -----------------------------
# Mean + SE labels (no group names)
# -----------------------------
make_stats_labels <- function(df, var) {
  df %>%
    group_by(group) %>%
    summarise(
      mean = mean(.data[[var]], na.rm = TRUE),
      se   = sd(.data[[var]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[var]]))),
      .groups = "drop"
    ) %>%
    mutate(label = sprintf("mean = %.2f (SE = %.2f)", mean, se)) %>%
    arrange(match(group, names(group_cols)))
}

stats_prev <- make_stats_labels(data, "prev")
stats_rho  <- make_stats_labels(data, "rho")

# -----------------------------
# Annotation position:
# right of mean lines + inside panel
# -----------------------------
annot_pos_right_of_mean_density <- function(df, var, means_df, bins = 35,
                                            y_start = 0.92, line_gap = 0.10,
                                            offset_frac = 0.06) {
  v <- df[[var]]
  x_max <- max(v, na.rm = TRUE)
  x_min <- min(v, na.rm = TRUE)
  x_rng <- x_max - x_min

  # x: just to the right of the largest mean line
  x <- max(means_df$mu, na.rm = TRUE) + offset_frac * x_rng
  x <- min(x, x_max - 0.02 * x_rng)  # keep inside panel

  # y: based on max density peak across groups
  y_top <- max(sapply(split(v, df$group), function(z) {
    h <- hist(z, breaks = bins, plot = FALSE)
    max(h$density)
  }), na.rm = TRUE)

  n_lines <- length(unique(df$group)) + 1  # +1 for p line
  y <- y_top * (y_start - (0:(n_lines - 1)) * line_gap)

  list(x = x, y = y)
}

pos_prev <- annot_pos_right_of_mean_density(data, "prev", means_prev, bins = 35)
pos_rho  <- annot_pos_right_of_mean_density(data, "rho",  means_rho,  bins = 35)

# -----------------------------
# Shift the whole block right so mean/SE aligns with p-line
# -----------------------------
block_shift_frac <- 0.10  # increase/decrease to taste

prev_rng <- diff(range(data$prev, na.rm = TRUE))
rho_rng  <- diff(range(data$rho,  na.rm = TRUE))

pos_prev$x <- min(pos_prev$x + block_shift_frac * prev_rng,
                  max(data$prev, na.rm = TRUE) - 0.02 * prev_rng)

pos_rho$x  <- min(pos_rho$x + block_shift_frac * rho_rng,
                  max(data$rho,  na.rm = TRUE) - 0.02 * rho_rng)

# -----------------------------
# Build annotation data (ALL lines share the same x)
# -----------------------------
ann_prev <- stats_prev %>% mutate(x = pos_prev$x, y = pos_prev$y[seq_len(n())])
ann_rho  <- stats_rho  %>% mutate(x = pos_rho$x,  y = pos_rho$y[seq_len(n())])

p_prev_df <- data.frame(
  x = pos_prev$x,
  y = pos_prev$y[nrow(stats_prev) + 1],
  label = fmt_p(p_prev)
)

p_rho_df <- data.frame(
  x = pos_rho$x,
  y = pos_rho$y[nrow(stats_rho) + 1],
  label = fmt_p(p_rho)
)

# -----------------------------
# Plot 1: prev (density)
# -----------------------------
p1 <- ggplot(data, aes(x = prev, fill = group)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 35, alpha = 0.5, color = "white", position = "identity") +
  geom_vline(data = means_prev, aes(xintercept = mu, color = group),
             linewidth = 1.1, show.legend = FALSE) +
  geom_text(
    data = ann_prev,
    aes(x = x, y = y, label = label, color = group),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4,
    fontface = "bold",
    show.legend = FALSE
  ) +
  geom_text(
    data = p_prev_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4,
	fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols) +
  labs(x = "Prevalence Percentage Drop", y = "Density") +
  theme_minimal()

# -----------------------------
# Plot 2: rho (density)
# -----------------------------
p2 <- ggplot(data, aes(x = rho, fill = group)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 35, alpha = 0.5, color = "white", position = "identity") +
  geom_vline(data = means_rho, aes(xintercept = mu, color = group),
             linewidth = 1.1, show.legend = FALSE) +
  geom_text(
    data = ann_rho,
    aes(x = x, y = y, label = label, color = group),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4,
    fontface = "bold",
    show.legend = FALSE
  ) +
  geom_text(
    data = p_rho_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4,
	fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols) +
  labs(x = "Heritability Percentage Drop", y = "") +
  theme_minimal()

# -----------------------------
# Combine + legend fix
# -----------------------------
combined_plot <- (p1 + p2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "",
    subtitle = paste0(
      "\n",
      "Settings: Prevalence = ", prev, "%, Heritability = ", h2, "%"
    )
  ) &
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  )

combined_plot <- combined_plot & labs(fill = NULL)

combined_plot

# -----------------------------
# Save (taller so nothing gets squeezed/clipped)
# -----------------------------
ggsave(
  filename = paste0("CIF_simulation_prev_", prev, "_h2_", h2, ".pdf"),
  plot = combined_plot,
  width = 12,
  height = 5,
  units = "in",
  limitsize = FALSE
)

ggsave(
  filename = paste0("CIF_simulation_prev_", prev, "_h2_", h2, ".png"),
  plot = combined_plot,
  width = 12,
  height = 5,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

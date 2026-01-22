#!/usr/bin/env Rscript
# =============================================================================
# BASELINE SENSITIVITY VISUALIZATION
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# Paths
output_dir <- "/data2/GWProject/outputs/Sensitivity/Baseline"

# =============================================================================
# Load Data
# =============================================================================

summary_df <- read.csv(file.path(output_dir, "baseline_sensitivity_summary.csv"))
aquifer_df <- read.csv(file.path(output_dir, "baseline_sensitivity_aquifer_comparison.csv"))

# =============================================================================
# Plot 1: Total Storage by Baseline Period
# =============================================================================

p1 <- ggplot(summary_df, aes(x = reorder(Label, Total_Storage_km3), y = Total_Storage_km3)) +
  geom_col(fill = "#2166AC", alpha = 0.8, width = 0.7) +
  geom_errorbar(aes(ymin = Total_Storage_km3 - 1.96*Total_SD_km3,
                    ymax = Total_Storage_km3 + 1.96*Total_SD_km3),
                width = 0.2, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  coord_flip() +
  labs(title = "Baseline Period Sensitivity Analysis",
       subtitle = "Total CONUS groundwater storage change (1940-2025 vs baseline)",
       x = NULL,
       y = expression("Storage Change (km"^3*")"),
       caption = "Error bars: 95% CI from Monte Carlo uncertainty") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "01_baseline_sensitivity_total.pdf"),
       p1, width = 8, height = 5, dpi = 300)

# =============================================================================
# Plot 2: Top 10 Aquifer Sensitivity
# =============================================================================

top10 <- aquifer_df %>%
  arrange(desc(Range_km3)) %>%
  head(10) %>%
  select(Aquifer, starts_with("Storage")) %>%
  pivot_longer(cols = starts_with("Storage"),
               names_to = "Period", values_to = "Storage_km3") %>%
  mutate(Period = gsub("Storage_", "", Period),
         Period = gsub("_", "-", Period))

p2 <- ggplot(top10, aes(x = Period, y = Storage_km3, group = Aquifer, color = Aquifer)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  scale_color_viridis_d(option = "turbo") +
  labs(title = "Aquifer Sensitivity to Baseline Period",
       subtitle = "Top 10 most sensitive aquifers",
       x = "Baseline Period",
       y = expression("Storage Change (km"^3*")"),
       color = "Aquifer") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40"),
    legend.position = "right",
    legend.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(output_dir, "02_aquifer_sensitivity.pdf"),
       p2, width = 10, height = 6, dpi = 300)

# =============================================================================
# Plot 3: Sensitivity Range Bar Chart
# =============================================================================

top_aquifers <- aquifer_df %>%
  arrange(desc(Range_km3)) %>%
  head(15) %>%
  mutate(Aquifer = gsub("_", " ", Aquifer))

p3 <- ggplot(top_aquifers, aes(x = reorder(Aquifer, Range_km3), y = Range_km3)) +
  geom_col(fill = "#D6604D", alpha = 0.8, width = 0.7) +
  coord_flip() +
  labs(title = "Aquifer Sensitivity to Baseline Period Choice",
       subtitle = "Range of storage estimates across 5 baseline periods",
       x = NULL,
       y = expression("Storage Range (km"^3*")")) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40"),
    panel.grid.major.y = element_blank()
  )

ggsave(file.path(output_dir, "03_sensitivity_range.pdf"),
       p3, width = 8, height = 6, dpi = 300)

# =============================================================================
# Summary Table
# =============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse=""), "\n")
cat("BASELINE SENSITIVITY SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse=""), "\n")

cat("\nTotal CONUS Storage by Baseline Period:\n")
cat("-" |> rep(50) |> paste(collapse=""), "\n")
for (i in 1:nrow(summary_df)) {
  row <- summary_df[i,]
  cat(sprintf("  %-25s %7.0f ± %4.0f km³\n",
              row$Label, row$Total_Storage_km3, row$Total_SD_km3))
}

# Key metrics
range_km3 <- max(summary_df$Total_Storage_km3) - min(summary_df$Total_Storage_km3)
mean_storage <- mean(summary_df$Total_Storage_km3)
range_pct <- 100 * range_km3 / abs(mean_storage)

cat("\n")
cat("-" |> rep(50) |> paste(collapse=""), "\n")
cat("KEY FINDINGS:\n")
cat("-" |> rep(50) |> paste(collapse=""), "\n")
cat(sprintf("  Storage range:       %.0f km³\n", range_km3))
cat(sprintf("  Range as %% of mean:  %.1f%%\n", range_pct))
cat(sprintf("  Mean across periods: %.0f km³\n", mean_storage))
cat(sprintf("  SD across periods:   %.0f km³\n", sd(summary_df$Total_Storage_km3)))

cat("\n")
cat("-" |> rep(50) |> paste(collapse=""), "\n")
cat("INTERPRETATION:\n")
cat("-" |> rep(50) |> paste(collapse=""), "\n")
cat("  Earlier baselines (1930-1945) show LARGER depletion\n")
cat("  because they capture more of pre-development conditions.\n")
cat("\n")
cat("  Later baselines (1950-1965) show SMALLER depletion\n")
cat("  because significant pumping had already begun.\n")
cat("\n")
cat("  The 1940-1955 baseline is a reasonable compromise:\n")
cat("  - Good data coverage across most aquifers\n")
cat("  - Pre-major irrigation expansion (High Plains)\n")
cat("  - Consistent with USGS historical assessments\n")

cat("\n")
cat("Plots saved to:\n")
cat(sprintf("  %s/01_baseline_sensitivity_total.pdf\n", output_dir))
cat(sprintf("  %s/02_aquifer_sensitivity.pdf\n", output_dir))
cat(sprintf("  %s/03_sensitivity_range.pdf\n", output_dir))

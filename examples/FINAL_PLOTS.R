#!/usr/bin/env Rscript
# =============================================================================
# FINAL ANALYSIS PLOTS
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(sf)
  library(cowplot)
  library(scales)
})

# Paths
results_file <- "/data2/GWProject/outputs/Phase3_Storage/STORAGE_RESULTS.csv"
sensitivity_file <- "/data2/GWProject/outputs/Phase4_Sensitivity/baseline_sensitivity.csv"
trends_file <- "/data2/GWProject/outputs/Phase5_Trends/trend_breakpoints.csv"
aquifer_gpkg <- "/data2/GWProject/data/EPA_Aquifers_Oct2025.gpkg"
output_dir <- "/data2/GWProject/outputs/Final/Plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
cat("Loading data...\n")
results <- read.csv(results_file)
sensitivity <- read.csv(sensitivity_file)
trends <- read.csv(trends_file)

# US states
us_states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))

# =============================================================================
# PLOT 1: Storage Change Bar Chart (Horizontal)
# =============================================================================

cat("Creating storage bar chart...\n")

results_plot <- results %>%
  mutate(Aquifer_Label = gsub("_", " ", Aquifer)) %>%
  arrange(Storage_km3)

p1 <- ggplot(results_plot, aes(x = reorder(Aquifer_Label, Storage_km3), y = Storage_km3)) +
  geom_col(aes(fill = Storage_km3 < 0), width = 0.75) +
  geom_errorbar(aes(ymin = Storage_km3 - 1.96*Storage_SD_km3,
                    ymax = Storage_km3 + 1.96*Storage_SD_km3),
                width = 0.3, linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray30") +
  scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"), guide = "none") +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Groundwater Storage Change by Aquifer",
    subtitle = "1940-2025 relative to 1940-1955 baseline",
    x = NULL,
    y = expression("Storage Change (km"^3*")"),
    caption = "Error bars: 95% CI | Negative (red) = depletion"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "01_storage_by_aquifer.pdf"), p1, width = 10, height = 9, dpi = 300)
ggsave(file.path(output_dir, "01_storage_by_aquifer.png"), p1, width = 10, height = 9, dpi = 150)

# =============================================================================
# PLOT 2: Spatial Map of Storage Change
# =============================================================================

cat("Creating spatial map...\n")

# Load aquifer boundaries
aquifers_sf <- st_read(aquifer_gpkg, quiet = TRUE)
aquifers_sf <- aquifers_sf %>%
  mutate(Aquifer_ID = gsub("[^A-Za-z0-9]", "_", AQ_NAME)) %>%
  mutate(Aquifer_ID = gsub("__+", "_", Aquifer_ID))

# Join with results
aquifers_data <- aquifers_sf %>%
  left_join(results, by = c("Aquifer_ID" = "Aquifer"))

p2 <- ggplot() +
  geom_sf(data = us_states, fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_sf(data = aquifers_data %>% filter(!is.na(Storage_km3)),
          aes(fill = Storage_km3), color = "gray30", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC",
    midpoint = 0,
    name = expression("km"^3),
    limits = c(-700, 100),
    oob = squish
  ) +
  coord_sf(xlim = c(-125, -66), ylim = c(24, 50), crs = 4326) +
  labs(
    title = "Groundwater Storage Change (1940-2025)",
    subtitle = "Relative to 1940-1955 baseline",
    caption = "Red = depletion | Blue = recovery"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40"),
    legend.position = "right",
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )

ggsave(file.path(output_dir, "02_storage_map.pdf"), p2, width = 12, height = 8, dpi = 300)
ggsave(file.path(output_dir, "02_storage_map.png"), p2, width = 12, height = 8, dpi = 150)

# =============================================================================
# PLOT 3: Baseline Sensitivity
# =============================================================================

cat("Creating sensitivity plot...\n")

p3 <- ggplot(sensitivity, aes(x = Baseline_Period, y = Total_Storage_km3)) +
  geom_col(fill = "#4575B4", alpha = 0.8, width = 0.7) +
  geom_errorbar(aes(ymin = Total_Storage_km3 - 1.96*Total_SD_km3,
                    ymax = Total_Storage_km3 + 1.96*Total_SD_km3),
                width = 0.2, linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Sensitivity to Baseline Period",
    subtitle = "Total CONUS storage change for different baselines",
    x = "Baseline Period",
    y = expression("Total Storage Change (km"^3*")"),
    caption = "Error bars: 95% CI"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(output_dir, "03_baseline_sensitivity.pdf"), p3, width = 8, height = 6, dpi = 300)

# =============================================================================
# PLOT 4: Trend Patterns Map
# =============================================================================

cat("Creating trend pattern map...\n")

# Join trends with aquifers
aquifers_trends <- aquifers_sf %>%
  left_join(trends, by = c("Aquifer_ID" = "Aquifer"))

pattern_colors <- c(
  "ACCELERATING" = "#D73027",
  "DECELERATING" = "#4575B4",
  "RECOVERY" = "#1A9850",
  "STABLE" = "#FFFFBF",
  "LINEAR" = "#FEE090"
)

p4 <- ggplot() +
  geom_sf(data = us_states, fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_sf(data = aquifers_trends %>% filter(!is.na(Pattern)),
          aes(fill = Pattern), color = "gray30", linewidth = 0.3) +
  scale_fill_manual(values = pattern_colors, name = "Pattern") +
  coord_sf(xlim = c(-125, -66), ylim = c(24, 50), crs = 4326) +
  labs(
    title = "Groundwater Trend Patterns",
    subtitle = "Based on breakpoint analysis",
    caption = "Accelerating = depletion rate increasing"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40"),
    legend.position = "right",
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )

ggsave(file.path(output_dir, "04_trend_patterns_map.pdf"), p4, width = 12, height = 8, dpi = 300)

# =============================================================================
# PLOT 5: Top 10 Time Series
# =============================================================================

cat("Creating time series plots...\n")

# Load annual data for top depleted aquifers
top_aquifers <- results %>%
  arrange(Storage_km3) %>%
  head(6) %>%
  pull(Aquifer)

annual_dir <- "/data2/GWProject/outputs/Phase3_Storage/annual"
annual_data <- data.frame()

for (aq in top_aquifers) {
  f <- file.path(annual_dir, paste0(aq, "_annual.csv"))
  if (file.exists(f)) {
    df <- read.csv(f)
    df$Aquifer <- gsub("_", " ", aq)
    annual_data <- rbind(annual_data, df)
  }
}

p5 <- ggplot(annual_data, aes(x = Year, y = IDW_Mean, color = Aquifer)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 0.8, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_viridis_d(option = "turbo") +
  labs(
    title = "Water Level Anomaly Time Series",
    subtitle = "Top 6 most depleted aquifers (IDW-corrected)",
    x = "Year",
    y = "Anomaly (m)",
    caption = "Negative = deeper water table = depletion"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40"),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2))

ggsave(file.path(output_dir, "05_timeseries_top6.pdf"), p5, width = 10, height = 6, dpi = 300)

# =============================================================================
# PLOT 6: Summary Dashboard
# =============================================================================

cat("Creating dashboard...\n")

# Summary stats for annotation
total_storage <- sum(results$Storage_km3)
total_sd <- sqrt(sum(results$Storage_SD_km3^2))
n_depleted <- sum(results$Storage_km3 < 0)
n_aquifers <- nrow(results)

summary_text <- paste0(
  "Total CONUS: ", format(round(total_storage), big.mark = ","), " ± ",
  round(total_sd), " km³\n",
  n_depleted, "/", n_aquifers, " aquifers depleted"
)

# Create smaller versions for dashboard
p1_small <- p1 + theme(legend.position = "none", plot.title = element_text(size = 11))
p2_small <- p2 + theme(legend.key.height = unit(0.4, "cm"), plot.title = element_text(size = 11))
p3_small <- p3 + theme(plot.title = element_text(size = 11))
p5_small <- p5 + theme(legend.position = "right", plot.title = element_text(size = 11))

dashboard <- plot_grid(
  p2_small, p1_small,
  p3_small, p5_small,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D"),
  rel_widths = c(1.2, 1)
)

# Add title
title <- ggdraw() +
  draw_label("CONUS Groundwater Storage Analysis (1940-2025)",
             fontface = "bold", size = 16, x = 0.5) +
  draw_label(summary_text, size = 10, x = 0.5, y = 0.2, color = "gray30")

dashboard_final <- plot_grid(title, dashboard, ncol = 1, rel_heights = c(0.1, 1))

ggsave(file.path(output_dir, "06_dashboard.pdf"), dashboard_final, width = 14, height = 12, dpi = 300)
ggsave(file.path(output_dir, "06_dashboard.png"), dashboard_final, width = 14, height = 12, dpi = 150)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat(rep("=", 60), "\n", sep = "")
cat("FINAL PLOTS CREATED\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nTotal CONUS Storage Change:", round(total_storage), "±", round(total_sd), "km³\n")
cat("Aquifers analyzed:", n_aquifers, "\n")
cat("Depleted:", n_depleted, "\n")
cat("\nPlots saved to:", output_dir, "\n")
cat("  01_storage_by_aquifer.pdf/png\n")
cat("  02_storage_map.pdf/png\n")
cat("  03_baseline_sensitivity.pdf\n")
cat("  04_trend_patterns_map.pdf\n")
cat("  05_timeseries_top6.pdf\n")
cat("  06_dashboard.pdf/png\n")

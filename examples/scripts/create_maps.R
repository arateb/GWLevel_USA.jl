#!/usr/bin/env Rscript
# =============================================================================
# SPATIAL MAPS FOR GROUNDWATER ANALYSIS
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(dplyr)
  library(viridis)
  library(cowplot)
})

# Paths
aquifer_gpkg <- "/data2/GWProject/data/EPA_Aquifers_Oct2025.gpkg"
results_file <- "/data2/GWProject/outputs/BayesianSpatial/BAYESIAN_STORAGE_RESULTS.csv"
breakpoints_file <- "/data2/GWProject/outputs/TrendAnalysis/trend_breakpoints.csv"
output_dir <- "/data2/GWProject/outputs/Maps"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Load Data
# =============================================================================

cat("Loading data...\n")

# Load aquifer boundaries
aquifers <- st_read(aquifer_gpkg, quiet = TRUE)

# Standardize aquifer names for joining
aquifers <- aquifers %>%
  mutate(Aquifer_ID = gsub("[^A-Za-z0-9]", "_", AQ_NAME))

# Load results
results <- read.csv(results_file)
breakpoints <- read.csv(breakpoints_file)

# Join results to aquifers
aquifers_data <- aquifers %>%
  left_join(results, by = c("Aquifer_ID" = "Aquifer")) %>%
  left_join(breakpoints %>% select(Aquifer, Pattern, Breakpoint, Acceleration_m_dec2),
            by = c("Aquifer_ID" = "Aquifer"))

# US states for background
us_states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))

cat("Data loaded successfully\n")

# =============================================================================
# Map 1: Storage Change by Aquifer
# =============================================================================

cat("Creating storage change map...\n")

p1 <- ggplot() +
  geom_sf(data = us_states, fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_sf(data = aquifers_data %>% filter(!is.na(Storage_km3)),
          aes(fill = Storage_km3), color = "gray30", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC",
    midpoint = 0,
    name = expression("Storage (km"^3*")"),
    limits = c(-700, 100),
    oob = scales::squish
  ) +
  coord_sf(xlim = c(-125, -66), ylim = c(24, 50), crs = 4326) +
  labs(
    title = "Groundwater Storage Change (1940-2025)",
    subtitle = "Relative to 1940-1955 baseline",
    caption = "Negative (red) = depletion | Positive (blue) = recovery"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40"),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    panel.grid = element_blank()
  )

ggsave(file.path(output_dir, "01_storage_change_map.pdf"),
       p1, width = 12, height = 8, dpi = 300)
ggsave(file.path(output_dir, "01_storage_change_map.png"),
       p1, width = 12, height = 8, dpi = 150)

# =============================================================================
# Map 2: Depletion Rate (m/decade)
# =============================================================================

cat("Creating depletion rate map...\n")

p2 <- ggplot() +
  geom_sf(data = us_states, fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_sf(data = aquifers_data %>% filter(!is.na(Trend_m_decade)),
          aes(fill = Trend_m_decade), color = "gray30", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC",
    midpoint = 0,
    name = "Rate (m/decade)",
    limits = c(-3, 1),
    oob = scales::squish
  ) +
  coord_sf(xlim = c(-125, -66), ylim = c(24, 50), crs = 4326) +
  labs(
    title = "Groundwater Depletion Rate",
    subtitle = "Linear trend in water level anomaly",
    caption = "Negative (red) = declining water table"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40"),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    panel.grid = element_blank()
  )

ggsave(file.path(output_dir, "02_depletion_rate_map.pdf"),
       p2, width = 12, height = 8, dpi = 300)

# =============================================================================
# Map 3: Trend Pattern (Accelerating vs Decelerating)
# =============================================================================

cat("Creating trend pattern map...\n")

# Simplify patterns
aquifers_data <- aquifers_data %>%
  mutate(Pattern_Simple = case_when(
    Pattern == "ACCELERATING_DEPLETION" ~ "Accelerating",
    Pattern == "DECELERATING_DEPLETION" ~ "Decelerating",
    Pattern == "RECOVERY" ~ "Recovery",
    Pattern == "REVERSAL_TO_DEPLETION" ~ "Reversal",
    TRUE ~ "Stable/Linear"
  ))

pattern_colors <- c(
  "Accelerating" = "#D73027",
  "Decelerating" = "#4575B4",
  "Recovery" = "#1A9850",
  "Reversal" = "#F46D43",
  "Stable/Linear" = "#FFFFBF"
)

p3 <- ggplot() +
  geom_sf(data = us_states, fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_sf(data = aquifers_data %>% filter(!is.na(Pattern_Simple)),
          aes(fill = Pattern_Simple), color = "gray30", linewidth = 0.3) +
  scale_fill_manual(values = pattern_colors, name = "Trend Pattern") +
  coord_sf(xlim = c(-125, -66), ylim = c(24, 50), crs = 4326) +
  labs(
    title = "Groundwater Trend Patterns",
    subtitle = "Detected from breakpoint analysis",
    caption = "Accelerating = depletion rate increasing over time"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40"),
    legend.position = "bottom",
    panel.grid = element_blank()
  )

ggsave(file.path(output_dir, "03_trend_pattern_map.pdf"),
       p3, width = 12, height = 8, dpi = 300)

# =============================================================================
# Map 4: Uncertainty (Storage SD)
# =============================================================================

cat("Creating uncertainty map...\n")

p4 <- ggplot() +
  geom_sf(data = us_states, fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_sf(data = aquifers_data %>% filter(!is.na(Storage_SD_km3)),
          aes(fill = Storage_SD_km3), color = "gray30", linewidth = 0.3) +
  scale_fill_viridis_c(
    option = "plasma",
    name = expression("SD (km"^3*")"),
    trans = "sqrt"
  ) +
  coord_sf(xlim = c(-125, -66), ylim = c(24, 50), crs = 4326) +
  labs(
    title = "Storage Estimate Uncertainty",
    subtitle = "Monte Carlo standard deviation",
    caption = "Higher values indicate greater uncertainty"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40"),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    panel.grid = element_blank()
  )

ggsave(file.path(output_dir, "04_uncertainty_map.pdf"),
       p4, width = 12, height = 8, dpi = 300)

# =============================================================================
# Combined Dashboard
# =============================================================================

cat("Creating combined dashboard...\n")

dashboard <- plot_grid(
  p1 + theme(legend.position = "right"),
  p2 + theme(legend.position = "right"),
  p3 + theme(legend.position = "right"),
  p4 + theme(legend.position = "right"),
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D")
)

ggsave(file.path(output_dir, "05_map_dashboard.pdf"),
       dashboard, width = 16, height = 14, dpi = 300)

# =============================================================================
# Summary
# =============================================================================

cat("\n")
cat(rep("=", 60), "\n", sep = "")
cat("MAPS CREATED\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nOutput directory:", output_dir, "\n")
cat("\nFiles:\n")
cat("  01_storage_change_map.pdf/png\n")
cat("  02_depletion_rate_map.pdf\n")
cat("  03_trend_pattern_map.pdf\n")
cat("  04_uncertainty_map.pdf\n")
cat("  05_map_dashboard.pdf\n")

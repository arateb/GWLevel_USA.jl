# Tutorial: Getting Started with GWLevel_USA

This tutorial walks through a complete groundwater storage analysis using the GWLevel_USA package.

---

## Installation

```julia
using Pkg

# From local path
Pkg.develop(path="/data/files/pkgs/GWLevel_USA")

# Or from GitHub (when published)
# Pkg.add(url="https://github.com/arateb/GWLevel_USA.jl")
```

---

## Quick Start

```julia
using GWLevel_USA

# Load your data
data = load_harmonized_data("/path/to/groundwater_data.csv")

# Define aquifer areas (km²)
areas = Dict(
    "High_Plains" => 505486.0,
    "Central_Valley" => 62421.0
)

# Run analysis
results = run_storage_analysis(data; aquifer_areas=areas)

# View summary
print_summary(results)
```

---

## Step-by-Step Guide

### Step 1: Prepare Your Data

Your CSV file needs these columns:

| Column | Type | Description |
|--------|------|-------------|
| WellID | String | Unique well identifier |
| Latitude | Float | WGS84 latitude (degrees) |
| Longitude | Float | WGS84 longitude (degrees) |
| Year | Int | Measurement year |
| DepthToWater_m | Float | Depth to water (meters) |
| Source | String | Data source (e.g., "NWIS") |
| Aquifer | String | Aquifer name |
| Season | String | Season (optional) |

**Example data format:**
```csv
WellID,Latitude,Longitude,Year,DepthToWater_m,Source,Aquifer,Season
TX001,32.5,-101.2,1950,15.3,TX,High_Plains,Spring
TX001,32.5,-101.2,1960,18.7,TX,High_Plains,Spring
TX001,32.5,-101.2,1970,25.1,TX,High_Plains,Spring
```

### Step 2: Load and Filter Data

```julia
using GWLevel_USA

# Load data
data = load_harmonized_data("/data/groundwater.csv")

# Check what we loaded
println("Records: ", nrow(data))
println("Wells: ", length(unique(data.WellID)))
println("Aquifers: ", unique(data.Aquifer))

# Filter to pre-irrigation seasons (recommended)
data = filter_by_season(data, ["Winter", "Spring"])
```

### Step 3: Configure Analysis

```julia
# Use default config
cfg = GWConfig()

# Or customize
cfg = GWConfig(
    baseline_start = 1940,      # Baseline period start
    baseline_end = 1955,        # Baseline period end
    min_obs_per_well = 5,       # Minimum observations per well
    idw_power = 2.0,            # IDW power parameter
    mc_samples = 10000          # Monte Carlo samples
)
```

### Step 4: Define Aquifer Areas

Aquifer areas must be in km², computed using an equal-area projection (EPSG:5070 recommended).

```julia
aquifer_areas = Dict{String, Float64}(
    "High_Plains" => 505486.0,
    "C_S_High_Plains" => 203400.0,
    "N_High_Plains" => 248320.0,
    "Basin_and_Range_Basin_Fill" => 429865.0,
    "Central_Valley" => 62421.0,
    # ... add more aquifers
)
```

### Step 5: Run Analysis

```julia
# Run for all aquifers
results = run_storage_analysis(data; cfg=cfg, aquifer_areas=aquifer_areas)

# Or analyze single aquifer
hp_data = get_wells_by_aquifer(data, "High_Plains")
hp_result = run_aquifer_analysis(hp_data, "High_Plains", 505486.0; cfg=cfg)
```

### Step 6: Examine Results

```julia
# Print summary
print_summary(results)

# Get summary table
summary_df = summarize_results(results)
println(summary_df)

# Access individual result
for r in results
    println("$(r.aquifer): $(round(r.storage.storage_mean_km3, digits=1)) ± $(round(r.storage.storage_sd_km3, digits=1)) km³")
end
```

### Step 7: Export Results

```julia
# Export to CSV files
export_results(results, "/output/storage_results/")

# Files created:
# - storage_summary.csv
# - total_conus.csv
# - annual/<Aquifer>_annual.csv
```

---

## Working with Individual Components

### Compute Baselines Manually

```julia
# Load data for one aquifer
data = load_harmonized_data("data.csv")
hp_data = get_wells_by_aquifer(data, "High_Plains")

# Compute baselines
baselines = compute_baseline_all(hp_data; cfg=GWConfig())
println("Wells with baseline: ", nrow(baselines))

# Add anomalies
hp_data = leftjoin(hp_data, baselines, on=:WellID)
hp_data.Anomaly_m = compute_anomaly.(hp_data.DepthToWater_m, hp_data.Baseline_m)
```

### Apply Quality Control

```julia
# Run Bayesian QC
hp_clean = run_bayesian_qc(hp_data)

# Check outlier rate
outlier_rate = mean(hp_clean.is_outlier)
println("Outlier rate: $(round(outlier_rate * 100, digits=1))%")

# Filter to inliers
hp_inliers = hp_clean[.!hp_clean.is_outlier, :]
```

### Spatial Interpolation

```julia
# Project coordinates
xs, ys = project_to_albers(hp_data.Longitude, hp_data.Latitude)

# IDW interpolation
z_pred = idw_interpolate(xs, ys, hp_data.Anomaly_m, [x_center], [y_center])

# Or get IDW mean directly
idw_mean = idw_interpolate_single(xs, ys, hp_data.Anomaly_m)
```

### Custom Storage Estimation

```julia
# Get Sy prior
sy = get_sy_prior("High_Plains")
println("Sy: $(sy.mean) ± $(sy.sd)")

# Estimate storage
result = estimate_storage(-15.0, 505486.0, sy; n_samples=10000)

println("Storage: $(round(result.storage_mean_km3)) ± $(round(result.storage_sd_km3)) km³")
println("95% CI: [$(round(result.ci_95[1])), $(round(result.ci_95[2]))]")
```

---

## Advanced Usage

### Custom Sy Priors

```julia
# Define custom prior
my_sy = SyPrior(0.12, 0.03, "My Study (2025)", :unconfined)

# Use in estimation
result = estimate_storage(-10.0, 100000.0, my_sy)
```

### Sensitivity Analysis

```julia
# Test different Sy values
for sy_mult in [0.8, 1.0, 1.2]
    sy = SyPrior(0.15 * sy_mult, 0.04, "Test", :unconfined)
    result = estimate_storage(-15.0, 200000.0, sy)
    println("Sy=$(round(sy.mean, digits=3)): Storage=$(round(result.storage_mean_km3)) km³")
end
```

### Trend Analysis

```julia
# Get annual time series
annual = results[1].annual_data

# Compute trend
trend = compute_trend(annual.Year, annual.MeanAnomaly_m)
println("Trend: $(round(trend.slope_decade, digits=2)) m/decade (R²=$(round(trend.r2, digits=2)))")
```

---

## Common Issues

### "No data remain after filtering"
- Check that your baseline period has data
- Reduce `min_obs_per_well` or `min_years_baseline`

### "Insufficient wells with baseline"
- Expand baseline period
- Check aquifer name spelling

### Large uncertainties
- Normal for aquifers with high Sy variability
- Check data quality (run QC first)

---

## Tips

1. **Always run QC first** - Outliers can significantly bias results
2. **Use Winter/Spring data** - Pre-irrigation measurements are more representative
3. **Check baseline coverage** - Wells need data in 1940-1955 period
4. **Validate against literature** - Compare with USGS reports
5. **Report uncertainties** - Always include SD or CI with estimates

---

*Last updated: January 2026*

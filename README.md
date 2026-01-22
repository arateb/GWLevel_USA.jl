# GWLevel_USA

[![Julia](https://img.shields.io/badge/Julia-1.9+-blue.svg)](https://julialang.org/)
[![License](https://img.shields.io/badge/License-Research-green.svg)]()

Bayesian framework for estimating groundwater storage changes across the conterminous United States using in-situ well observations.

**Author:** Dr. Ashraf Rateb, Bureau of Economic Geology, University of Texas at Austin

---

## Features

- **Data Processing**: Load and harmonize groundwater data from multiple sources (USGS NWIS, CA DWR, TX TWDB)
- **Spatial Analysis**: IDW interpolation to correct monitoring bias, equal-area projection (EPSG:5070)
- **Bayesian Estimation**: Monte Carlo uncertainty quantification with literature-based Sy priors
- **Quality Control**: t-mixture outlier detection, well classification
- **Flexible Configuration**: Customizable baseline periods, spatial parameters, and analysis settings

---

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/arateb/GWLevel_USA.jl")
```

Or for local development:

```julia
using Pkg
Pkg.activate("/path/to/GWLevel_USA")
Pkg.instantiate()
```

---

## Quick Start

```julia
using GWLevel_USA

# Load data
data = load_harmonized_data("/data/GW_CONUS_Harmonized.csv")

# Configure analysis
cfg = GWConfig(
    baseline_start = 1940,
    baseline_end = 1955,
    valid_seasons = ["Winter", "Spring"],
    mc_samples = 10000
)

# Define aquifer areas (km²)
aquifer_areas = Dict(
    "C_S_High_Plains" => 203400.0,
    "Basin_and_Range_Basin_Fill" => 429865.0,
    # ... more aquifers
)

# Run analysis
results = run_storage_analysis(data; cfg=cfg, aquifer_areas=aquifer_areas)

# Print summary
print_summary(results)

# Export results
export_results(results, "outputs/")
```

---

## Methodology

### 1. Baseline Calculation

Wells with observations during the baseline period (default: 1940-1955) are used to establish pre-development water levels:

```
Baseline_i = mean(Depth_i,t) for t ∈ [1940, 1955]
```

### 2. Anomaly Calculation

Water level anomaly follows GRACE convention (negative = depletion):

```
Anomaly = Baseline - Current_Depth
```

### 3. Spatial Weighting

IDW interpolation corrects for monitoring well clustering:

```
z_IDW = Σ(w_i × z_i) / Σ(w_i)
where w_i = 1 / d_i^p
```

### 4. Storage Estimation

Storage change with Monte Carlo uncertainty:

```
ΔS (km³) = Sy × Area (km²) × Δh (m) / 1000
```

Sy sampled from truncated normal distribution based on literature values.

---

## Package Structure

```
GWLevel_USA/
├── src/
│   ├── GWLevel_USA.jl      # Main module
│   ├── Config.jl           # Configuration types
│   ├── DataIO.jl           # Data loading
│   ├── Spatial.jl          # Projection, IDW
│   ├── Bayesian.jl         # Storage estimation
│   ├── QualityControl.jl   # Outlier detection
│   └── Analysis.jl         # High-level API
├── test/
│   └── runtests.jl
├── examples/
│   ├── run_analysis.jl     # Main CONUS analysis
│   └── scripts/
│       ├── baseline_sensitivity.jl   # Baseline period sensitivity
│       ├── trend_breakpoints.jl      # Piecewise trend analysis
│       ├── plot_baseline_sensitivity.R
│       └── create_maps.R             # Spatial visualization
├── docs/
├── Project.toml
└── README.md
```

---

## Analysis Scripts

### Baseline Sensitivity Analysis
Tests storage estimates across different baseline periods (1930-1945, 1935-1950, 1940-1955, 1945-1960, 1950-1965).

```bash
julia --project=. examples/scripts/baseline_sensitivity.jl
```

### Trend Breakpoint Analysis
Detects acceleration/deceleration periods using piecewise linear regression.

```bash
julia --project=. examples/scripts/trend_breakpoints.jl
```

### Spatial Maps (R)
Creates choropleth maps of storage change, depletion rates, and trend patterns.

```bash
Rscript examples/scripts/create_maps.R
```

---

## Key Functions

### Data Loading
```julia
load_harmonized_data(filepath; cfg)  # Load CSV data
filter_by_season(df, seasons)         # Filter by season
compute_well_statistics(df)           # Per-well stats
```

### Spatial Methods
```julia
project_to_albers(lon, lat)           # WGS84 → EPSG:5070
idw_interpolate(x, y, z, xp, yp)     # IDW interpolation
compute_nyquist_spacing(x, y)         # Optimal grid spacing
```

### Bayesian Estimation
```julia
compute_baseline(df, well_id)         # Baseline depth
compute_anomaly(depth, baseline)      # Water level anomaly
estimate_storage(anomaly, area, sy)   # Storage with MC
run_monte_carlo(...)                  # Direct MC sampling
```

### Quality Control
```julia
detect_outliers(values; z_threshold)  # MAD-based detection
run_bayesian_qc(df)                   # Full QC pipeline
classify_well(residuals, years)       # Well classification
```

### Analysis
```julia
run_aquifer_analysis(df, aquifer, area)  # Single aquifer
run_storage_analysis(data; cfg, areas)   # All aquifers
summarize_results(results)               # Summary table
export_results(results, output_dir)      # Export to CSV
```

---

## Configuration Options

```julia
GWConfig(
    # Temporal
    baseline_start = 1940,
    baseline_end = 1955,
    valid_seasons = ["Winter", "Spring"],
    year_min = 1940,
    year_max = 2025,

    # Data quality
    min_obs_per_well = 5,
    min_years_baseline = 2,
    min_wells_per_aquifer = 5,
    outlier_zscore = 3.5,

    # Spatial
    idw_power = 2.0,
    idw_neighbors = 12,
    grid_resolution = 0.25,

    # Monte Carlo
    mc_samples = 10000
)
```

---

## Specific Yield Priors

Literature-based Sy values are included for 30 major US aquifers:

| Aquifer | Sy | SD | Reference |
|---------|----|----|-----------|
| High Plains | 0.15 | 0.04 | Gutentag et al. (1984) |
| Central Valley | 0.12 | 0.04 | Faunt (2009) |
| Basin and Range | 0.18 | 0.04 | Pool & Coes (1999) |
| Mississippi Valley | 0.22 | 0.05 | Clark & Hart (2009) |
| Floridan (confined) | 0.001 | 0.0005 | Miller (1986) |

See `src/Config.jl` for complete list with references.

---

## Dependencies

- CSV.jl
- DataFrames.jl
- Distributions.jl
- NearestNeighbors.jl
- StatsBase.jl

---

## Testing

```julia
using Pkg
Pkg.test("GWLevel_USA")
```

---

## Citation

If you use this package, please cite:

```bibtex
@software{gwlevel_usa,
  author = {Rateb, Ashraf},
  title = {GWLevel_USA: Bayesian Groundwater Storage Estimation},
  year = {2026},
  institution = {Bureau of Economic Geology, UT Austin}
}
```

---

## References

- Gutentag, E.D., et al. (1984). USGS PP 1400-B.
- Faunt, C.C. (2009). USGS PP 1766.
- McGuire, V.L. (2017). USGS SIR 2017-5040.
- Pool, D.R., & Coes, A.L. (1999). USGS WRIR 99-4197.
- Konikow, L.F. (2013). USGS SIR 2013-5079.

---

## License

Research use. Contact author for collaboration.

---

## Key Results

Total CONUS groundwater storage change: **-2,232 ± 353 km³** (1940-2025 vs 1940-1955 baseline)

| Top Depleted Aquifers | Storage (km³) |
|-----------------------|---------------|
| Central/South High Plains | -623 ± 176 |
| Basin and Range | -594 ± 202 |
| Arizona Alluvials | -253 ± 82 |
| SE Coastal Plain | -145 ± 70 |
| Mississippi River Valley | -111 ± 35 |

---

*Last updated: January 2026*

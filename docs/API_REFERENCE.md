# API Reference

## Module: GWLevel_USA

Complete API documentation for the GWLevel_USA Julia package.

---

## Configuration

### `GWConfig`

Configuration parameters for groundwater analysis.

```julia
GWConfig(;
    # Temporal parameters
    baseline_start::Int = 1940,
    baseline_end::Int = 1955,
    valid_seasons::Vector{String} = ["Winter", "Spring"],
    year_min::Int = 1940,
    year_max::Int = 2025,
    recent_year::Int = 2020,

    # Data quality thresholds
    min_obs_per_well::Int = 5,
    min_years_baseline::Int = 2,
    min_wells_per_aquifer::Int = 5,
    min_records_per_aquifer::Int = 100,

    # Outlier detection
    outlier_zscore::Float64 = 3.5,
    outlier_prior::Float64 = 0.01,

    # Spatial parameters
    idw_power::Float64 = 2.0,
    idw_neighbors::Int = 12,
    grid_resolution::Float64 = 0.25,

    # Monte Carlo
    mc_samples::Int = 10000,

    # Projection (EPSG:5070)
    proj_lat0::Float64 = 23.0,
    proj_lon0::Float64 = -96.0,
    proj_lat1::Float64 = 29.5,
    proj_lat2::Float64 = 45.5
)
```

**Example:**
```julia
cfg = GWConfig(baseline_start=1950, baseline_end=1960, mc_samples=5000)
```

---

### `SyPrior`

Specific yield prior distribution.

```julia
SyPrior(mean::Float64, sd::Float64, reference::String, confinement::Symbol)
```

**Fields:**
- `mean`: Prior mean value
- `sd`: Prior standard deviation
- `reference`: Literature citation
- `confinement`: `:unconfined`, `:confined`, or `:mixed`

**Example:**
```julia
sy = SyPrior(0.15, 0.04, "Gutentag et al. (1984)", :unconfined)
```

---

### `get_sy_prior`

Get specific yield prior for an aquifer.

```julia
get_sy_prior(aquifer_name::String) -> SyPrior
```

Returns default `SyPrior(0.10, 0.05, "Default", :mixed)` if aquifer not found.

**Example:**
```julia
sy = get_sy_prior("C_S_High_Plains")
# SyPrior(0.15, 0.04, "Gutentag et al. (1984)", :unconfined)
```

---

## Data I/O

### `load_harmonized_data`

Load harmonized groundwater data from CSV.

```julia
load_harmonized_data(filepath::String; cfg::GWConfig=DEFAULT_CONFIG) -> DataFrame
```

**Required CSV columns:**
- `WellID`: Unique well identifier
- `Latitude`, `Longitude`: WGS84 coordinates
- `Year`: Measurement year
- `DepthToWater_m` (or `DepthToWater_ft`): Depth to water
- `Source`: Data source (NWIS, CA, TX)
- `Aquifer`: Aquifer name
- `Season` (optional): Season of measurement

**Example:**
```julia
data = load_harmonized_data("/data/GW_CONUS.csv")
```

---

### `filter_by_season`

Filter data to specified seasons.

```julia
filter_by_season(df::DataFrame, seasons::Vector{String}) -> DataFrame
```

**Example:**
```julia
winter_spring = filter_by_season(data, ["Winter", "Spring"])
```

---

### `filter_by_years`

Filter data to year range.

```julia
filter_by_years(df::DataFrame, year_start::Int, year_end::Int) -> DataFrame
```

---

### `compute_well_statistics`

Compute per-well summary statistics.

```julia
compute_well_statistics(df::DataFrame) -> DataFrame
```

**Returns DataFrame with:**
- `WellID`, `Latitude`, `Longitude`, `Aquifer`, `Source`
- `n_obs`, `n_years`, `year_min`, `year_max`
- `depth_mean`, `depth_median`, `depth_sd`, `depth_iqr`

---

### `get_wells_by_aquifer`

Extract data for a specific aquifer.

```julia
get_wells_by_aquifer(df::DataFrame, aquifer::String) -> DataFrame
```

---

### `get_aquifer_list`

Get list of unique aquifers.

```julia
get_aquifer_list(df::DataFrame) -> Vector{String}
```

---

## Spatial Methods

### `project_to_albers`

Project WGS84 to NAD83/Conus Albers (EPSG:5070).

```julia
# Single point
project_to_albers(lon::Real, lat::Real; cfg::GWConfig=DEFAULT_CONFIG) -> Tuple{Float64, Float64}

# Vectorized
project_to_albers(lons::AbstractVector, lats::AbstractVector; cfg::GWConfig=DEFAULT_CONFIG) -> Tuple{Vector{Float64}, Vector{Float64}}
```

**Returns:** `(x, y)` in meters

**Example:**
```julia
x, y = project_to_albers(-97.7, 30.3)
xs, ys = project_to_albers(df.Longitude, df.Latitude)
```

---

### `inverse_project`

Inverse projection from Albers to WGS84.

```julia
inverse_project(x::Real, y::Real; cfg::GWConfig=DEFAULT_CONFIG) -> Tuple{Float64, Float64}
```

**Returns:** `(lon, lat)` in degrees

---

### `idw_interpolate`

Inverse Distance Weighting interpolation.

```julia
idw_interpolate(
    x_obs::AbstractVector, y_obs::AbstractVector, z_obs::AbstractVector,
    x_pred::AbstractVector, y_pred::AbstractVector;
    idp::Real=2, nmax::Int=12, min_dist::Real=1.0
) -> Vector{Float64}
```

**Arguments:**
- `x_obs, y_obs`: Observation coordinates (meters)
- `z_obs`: Observation values
- `x_pred, y_pred`: Prediction coordinates
- `idp`: Power parameter (default: 2)
- `nmax`: Maximum neighbors (default: 12)

**Example:**
```julia
z_pred = idw_interpolate(x_wells, y_wells, anomaly, x_grid, y_grid; idp=2)
```

---

### `idw_interpolate_single`

Compute IDW mean at centroid.

```julia
idw_interpolate_single(x_obs, y_obs, z_obs; idp=2, nmax=12) -> Float64
```

---

### `compute_nyquist_spacing`

Compute optimal grid spacing from point distribution.

```julia
compute_nyquist_spacing(x::AbstractVector, y::AbstractVector) -> Float64
```

**Returns:** `median(NN_distance) / 2`

---

## Bayesian Estimation

### `compute_baseline`

Compute baseline depth for a well.

```julia
compute_baseline(df::DataFrame, well_id::Any;
                 baseline_start::Int=1940, baseline_end::Int=1955,
                 min_years::Int=2) -> Float64
```

**Returns:** Mean depth during baseline, or `NaN` if insufficient data.

---

### `compute_baseline_all`

Compute baselines for all wells.

```julia
compute_baseline_all(df::DataFrame; cfg::GWConfig=DEFAULT_CONFIG) -> DataFrame
```

**Returns:** DataFrame with `WellID` and `Baseline_m`

---

### `compute_anomaly`

Compute water level anomaly.

```julia
compute_anomaly(depth::Real, baseline::Real) -> Float64
```

**Convention:** `Anomaly = Baseline - Current` (negative = depletion)

---

### `estimate_storage`

Estimate storage with Monte Carlo uncertainty.

```julia
estimate_storage(
    anomaly_m::Real, area_km2::Real, sy::SyPrior;
    aquifer_name::String="Unknown",
    anomaly_sd_m::Real=0.0,
    n_samples::Int=10000
) -> StorageResult
```

**Returns:** `StorageResult` struct with mean, SD, CI, samples.

**Example:**
```julia
sy = get_sy_prior("High_Plains")
result = estimate_storage(-10.0, 500000.0, sy; n_samples=10000)
println("Storage: $(result.storage_mean_km3) ± $(result.storage_sd_km3) km³")
```

---

### `run_monte_carlo`

Simple Monte Carlo storage estimation.

```julia
run_monte_carlo(anomaly_m, anomaly_sd, sy_mean, sy_sd, area_km2;
                n_samples=10000) -> NamedTuple
```

**Returns:** `(mean, sd, q05, q25, q50, q75, q95)`

---

### `compute_trend`

Compute linear trend.

```julia
compute_trend(years::AbstractVector, values::AbstractVector) -> NamedTuple
```

**Returns:** `(slope_decade, intercept, r2)`

---

### `classify_quality`

Classify data quality tier.

```julia
classify_quality(n_years::Int, n_wells::Int, r2::Real, cv::Real) -> Symbol
```

**Returns:** `:HIGH`, `:MEDIUM`, or `:LOW`

---

## Quality Control

### `detect_outliers`

MAD-based outlier detection.

```julia
detect_outliers(values::AbstractVector; z_threshold::Real=3.5) -> BitVector
```

**Returns:** `true` where value is outlier.

---

### `robust_mad`

Median Absolute Deviation.

```julia
robust_mad(x::AbstractVector; constant::Real=1.4826) -> Float64
```

---

### `t_mixture_posterior`

Bayesian outlier probability.

```julia
t_mixture_posterior(residual, sigma, pi0;
                    df_in=7.0, df_out=1.0, scale_mult=3.0) -> Float64
```

**Returns:** P(outlier | residual)

---

### `run_bayesian_qc`

Full Bayesian QC pipeline.

```julia
run_bayesian_qc(df::DataFrame; cfg::GWConfig=DEFAULT_CONFIG) -> DataFrame
```

**Adds columns:** `fitted`, `residual`, `p_outlier`, `is_outlier`, `w_qc`

---

### `classify_well`

Classify well behavior.

```julia
classify_well(residuals, years; amp_small=1.0, amp_large=3.0,
              r2_min=0.70, n_min=3) -> WellClass
```

**Returns:** `LIKELY_CONFINED`, `LIKELY_UNCONFINED`, `INDETERMINATE`, etc.

---

## Analysis

### `run_aquifer_analysis`

Analyze single aquifer.

```julia
run_aquifer_analysis(df::DataFrame, aquifer::String, area_km2::Real;
                     cfg::GWConfig=DEFAULT_CONFIG) -> AquiferResult
```

---

### `run_storage_analysis`

Analyze all aquifers.

```julia
run_storage_analysis(data::DataFrame;
                     cfg::GWConfig=DEFAULT_CONFIG,
                     aquifer_areas::Dict{String,Float64}=Dict()) -> Vector{AquiferResult}
```

**Example:**
```julia
areas = Dict("High_Plains" => 505486.0, "Basin_and_Range" => 429865.0)
results = run_storage_analysis(data; cfg=cfg, aquifer_areas=areas)
```

---

### `summarize_results`

Create summary DataFrame.

```julia
summarize_results(results::Vector{AquiferResult}) -> DataFrame
```

---

### `export_results`

Export to CSV files.

```julia
export_results(results::Vector{AquiferResult}, output_dir::String)
```

**Creates:**
- `storage_summary.csv`
- `total_conus.csv`
- `annual/<Aquifer>_annual.csv`

---

### `print_summary`

Print console summary.

```julia
print_summary(results::Vector{AquiferResult})
```

---

## Result Types

### `StorageResult`

```julia
struct StorageResult
    aquifer::String
    area_km2::Float64
    anomaly_m::Float64
    anomaly_sd_m::Float64
    sy_prior::SyPrior
    sy_post_mean::Float64
    sy_post_sd::Float64
    storage_mean_km3::Float64
    storage_sd_km3::Float64
    storage_samples::Vector{Float64}
    ci_95::Tuple{Float64, Float64}
end
```

### `AquiferResult`

```julia
struct AquiferResult
    aquifer::String
    n_wells::Int
    n_years::Int
    year_range::Tuple{Int, Int}
    trend_m_decade::Float64
    trend_r2::Float64
    storage::StorageResult
    quality_tier::Symbol
    annual_data::DataFrame
end
```

---

*Last updated: January 2026*

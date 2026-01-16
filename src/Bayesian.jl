# =============================================================================
# Bayesian Storage Estimation
# =============================================================================

"""
    StorageResult

Results from Bayesian storage estimation for a single aquifer.

# Fields
- `aquifer::String`: Aquifer name
- `area_km2::Float64`: Aquifer area
- `anomaly_m::Float64`: Mean water level anomaly (m)
- `anomaly_sd_m::Float64`: Anomaly standard deviation
- `sy_prior::SyPrior`: Specific yield prior
- `sy_post_mean::Float64`: Posterior Sy mean
- `sy_post_sd::Float64`: Posterior Sy SD
- `storage_mean_km3::Float64`: Storage change mean (km³)
- `storage_sd_km3::Float64`: Storage change SD
- `storage_samples::Vector{Float64}`: MC samples (optional)
- `ci_95::Tuple{Float64, Float64}`: 95% credible interval
"""
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

"""
    AquiferResult

Complete analysis results for an aquifer.

# Fields
- `aquifer::String`: Aquifer name
- `n_wells::Int`: Number of wells
- `n_years::Int`: Number of years with data
- `year_range::Tuple{Int, Int}`: Year range
- `trend_m_decade::Float64`: Linear trend (m/decade)
- `trend_r2::Float64`: Trend R²
- `storage::StorageResult`: Storage estimation results
- `quality_tier::Symbol`: :HIGH, :MEDIUM, or :LOW
- `annual_data::DataFrame`: Annual time series
"""
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

"""
    compute_baseline(df::DataFrame, well_id::Any;
                     baseline_start::Int=1940, baseline_end::Int=1955) -> Float64

Compute baseline mean depth for a well.

# Arguments
- `df`: DataFrame with Year and DepthToWater_m columns
- `well_id`: Well identifier
- `baseline_start`, `baseline_end`: Baseline period

# Returns
Mean depth during baseline period, or NaN if insufficient data
"""
function compute_baseline(
    df::DataFrame,
    well_id::Any;
    baseline_start::Int=1940,
    baseline_end::Int=1955,
    min_years::Int=2
)
    # Filter to well and baseline period
    well_data = df[(df.WellID .== well_id) .&
                   (df.Year .>= baseline_start) .&
                   (df.Year .<= baseline_end), :]

    if nrow(well_data) < min_years
        return NaN
    end

    return mean(well_data.DepthToWater_m)
end

"""
    compute_baseline_all(df::DataFrame; cfg::GWConfig=DEFAULT_CONFIG) -> DataFrame

Compute baselines for all wells.

# Returns
DataFrame with WellID and Baseline_m columns
"""
function compute_baseline_all(df::DataFrame; cfg::GWConfig=DEFAULT_CONFIG)
    wells = unique(df.WellID)

    baselines = DataFrame(
        WellID = wells,
        Baseline_m = [compute_baseline(df, w;
                                       baseline_start=cfg.baseline_start,
                                       baseline_end=cfg.baseline_end,
                                       min_years=cfg.min_years_baseline) for w in wells]
    )

    # Remove wells without valid baseline
    baselines = baselines[.!isnan.(baselines.Baseline_m), :]

    @info "Computed baselines for $(nrow(baselines)) wells"
    return baselines
end

"""
    compute_anomaly(depth::Real, baseline::Real) -> Float64

Compute water level anomaly.

# Convention
Anomaly = Baseline - Current
- Negative anomaly = deeper water table = depletion (GRACE convention)
- Positive anomaly = shallower water table = recovery
"""
function compute_anomaly(depth::Real, baseline::Real)
    return baseline - depth
end

"""
    add_anomalies!(df::DataFrame, baselines::DataFrame)

Add anomaly column to data frame (modifies in place).
"""
function add_anomalies!(df::DataFrame, baselines::DataFrame)
    df = leftjoin(df, baselines, on=:WellID)
    df.Anomaly_m = compute_anomaly.(df.DepthToWater_m, df.Baseline_m)
    return df
end

"""
    estimate_storage(anomaly_m::Real, area_km2::Real, sy::SyPrior;
                     n_samples::Int=10000) -> StorageResult

Estimate storage change with Monte Carlo uncertainty.

# Arguments
- `anomaly_m`: Mean water level anomaly (m)
- `area_km2`: Aquifer area (km²)
- `sy`: Specific yield prior
- `n_samples`: Number of Monte Carlo samples

# Returns
StorageResult with mean, SD, and credible intervals

# Storage Equation
ΔS (km³) = Sy × Area (km²) × Δh (m) / 1000
"""
function estimate_storage(
    anomaly_m::Real,
    area_km2::Real,
    sy::SyPrior;
    aquifer_name::String="Unknown",
    anomaly_sd_m::Real=0.0,
    n_samples::Int=10000
)
    # Sample from truncated normal for Sy (bounded by physics)
    sy_dist = truncated(Normal(sy.mean, sy.sd), 0.0001, 0.5)
    sy_samples = rand(sy_dist, n_samples)

    # Sample anomaly if uncertainty provided
    if anomaly_sd_m > 0
        anomaly_samples = rand(Normal(anomaly_m, anomaly_sd_m), n_samples)
    else
        anomaly_samples = fill(anomaly_m, n_samples)
    end

    # Storage samples: ΔS = Sy × Area × Δh / 1000
    storage_samples = sy_samples .* area_km2 .* anomaly_samples ./ 1000.0

    # Summary statistics
    storage_mean = mean(storage_samples)
    storage_sd = std(storage_samples)
    ci_95 = (quantile(storage_samples, 0.025), quantile(storage_samples, 0.975))

    # Posterior Sy (just use prior for now - could update with data)
    sy_post_mean = mean(sy_samples)
    sy_post_sd = std(sy_samples)

    return StorageResult(
        aquifer_name,
        area_km2,
        anomaly_m,
        anomaly_sd_m,
        sy,
        sy_post_mean,
        sy_post_sd,
        storage_mean,
        storage_sd,
        storage_samples,
        ci_95
    )
end

"""
    run_monte_carlo(anomaly_m::Real, anomaly_sd::Real,
                    sy_mean::Real, sy_sd::Real,
                    area_km2::Real;
                    n_samples::Int=10000) -> NamedTuple

Simple Monte Carlo storage estimation.

# Returns
NamedTuple with mean, sd, q05, q25, q50, q75, q95
"""
function run_monte_carlo(
    anomaly_m::Real,
    anomaly_sd::Real,
    sy_mean::Real,
    sy_sd::Real,
    area_km2::Real;
    n_samples::Int=10000
)
    # Sample Sy (truncated normal)
    sy_dist = truncated(Normal(sy_mean, sy_sd), 0.0001, 0.5)
    sy_samples = rand(sy_dist, n_samples)

    # Sample anomaly
    if anomaly_sd > 0
        anomaly_samples = rand(Normal(anomaly_m, anomaly_sd), n_samples)
    else
        anomaly_samples = fill(anomaly_m, n_samples)
    end

    # Storage
    storage = sy_samples .* area_km2 .* anomaly_samples ./ 1000.0

    return (
        mean = mean(storage),
        sd = std(storage),
        q05 = quantile(storage, 0.05),
        q25 = quantile(storage, 0.25),
        q50 = quantile(storage, 0.50),
        q75 = quantile(storage, 0.75),
        q95 = quantile(storage, 0.95)
    )
end

"""
    compute_trend(years::AbstractVector, values::AbstractVector) -> NamedTuple

Compute linear trend using OLS.

# Returns
NamedTuple with slope (per decade), intercept, r2, pvalue
"""
function compute_trend(years::AbstractVector{<:Real}, values::AbstractVector{<:Real})
    n = length(years)
    if n < 3
        return (slope_decade=NaN, intercept=NaN, r2=NaN)
    end

    # Filter valid
    valid = isfinite.(years) .& isfinite.(values)
    t = years[valid]
    v = values[valid]
    n = length(t)

    if n < 3
        return (slope_decade=NaN, intercept=NaN, r2=NaN)
    end

    # OLS
    t_mean = mean(t)
    v_mean = mean(v)
    Stt = sum((t .- t_mean).^2)
    Stv = sum((t .- t_mean) .* (v .- v_mean))

    slope = Stt > 0 ? Stv / Stt : 0.0
    intercept = v_mean - slope * t_mean

    # R²
    fitted = intercept .+ slope .* t
    ss_res = sum((v .- fitted).^2)
    ss_tot = sum((v .- v_mean).^2)
    r2 = ss_tot > 0 ? 1.0 - ss_res / ss_tot : 0.0

    return (
        slope_decade = slope * 10,  # Convert to per decade
        intercept = intercept,
        r2 = r2
    )
end

"""
    classify_quality(n_years::Int, n_wells::Int, r2::Real, cv::Real) -> Symbol

Classify aquifer data quality tier.

# Criteria
- HIGH: CV < 0.5, R² > 0.5, Years ≥ 70, Wells ≥ 100
- MEDIUM: CV < 1.0, R² > 0.2, Years ≥ 50, Wells ≥ 30
- LOW: Everything else
"""
function classify_quality(n_years::Int, n_wells::Int, r2::Real, cv::Real)
    if cv < 0.5 && r2 > 0.5 && n_years >= 70 && n_wells >= 100
        return :HIGH
    elseif cv < 1.0 && r2 > 0.2 && n_years >= 50 && n_wells >= 30
        return :MEDIUM
    else
        return :LOW
    end
end

export compute_baseline_all, add_anomalies!, compute_trend, classify_quality

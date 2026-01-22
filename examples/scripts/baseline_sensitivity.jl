#!/usr/bin/env julia
# =============================================================================
# BASELINE PERIOD SENSITIVITY ANALYSIS
# =============================================================================
#
# Tests storage estimates across different baseline periods:
#   1. 1930-1945 (earliest available)
#   2. 1935-1950
#   3. 1940-1955 (default)
#   4. 1945-1960
#   5. 1950-1965
#
# =============================================================================

using Pkg
Pkg.activate("/data2/GWProject")

using DataFrames, CSV, Statistics, Dates, Random
using StatsBase: quantile, std, mad

Random.seed!(42)

include("/data2/GWProject/src/spatial_methods.jl")
using .SpatialMethods

# =============================================================================
# Configuration
# =============================================================================

const N_SAMPLES = 10000
const INPUT_DIR = "/data2/GWProject/outputs/ByAquifer"
const OUTPUT_DIR = "/data2/GWProject/outputs/Sensitivity/Baseline"
mkpath(OUTPUT_DIR)

# Baseline periods to test (start_year, end_year)
const BASELINE_PERIODS = [
    (1930, 1945, "Early (1930-1945)"),
    (1935, 1950, "Pre-war (1935-1950)"),
    (1940, 1955, "Default (1940-1955)"),
    (1945, 1960, "Post-war (1945-1960)"),
    (1950, 1965, "Late (1950-1965)")
]

# Sy priors (same as main analysis)
const SY_PRIORS = Dict(
    "N_High_Plains" => (sy_mean=0.15, sy_sd=0.04),
    "C_S_High_Plains" => (sy_mean=0.15, sy_sd=0.04),
    "High_Plains" => (sy_mean=0.15, sy_sd=0.04),
    "Central_Valley" => (sy_mean=0.12, sy_sd=0.04),
    "Basin_and_Range_Basin_Fill" => (sy_mean=0.18, sy_sd=0.04),
    "Arizona_Alluvials" => (sy_mean=0.18, sy_sd=0.05),
    "Mississippi_River_Valley_Alluvial" => (sy_mean=0.22, sy_sd=0.05),
    "Mississippi_Embayment" => (sy_mean=0.12, sy_sd=0.04),
    "Snake" => (sy_mean=0.08, sy_sd=0.04),
    "Columbia_Plateau_Basin_fill" => (sy_mean=0.12, sy_sd=0.05),
    "Columbia_Plateau_Basaltic_rock" => (sy_mean=0.05, sy_sd=0.03),
    "Coastal_Lowland" => (sy_mean=0.08, sy_sd=0.04),
    "Texas_Gulf_Coast" => (sy_mean=0.06, sy_sd=0.03),
    "N_Atlantic_Coastal_Plain" => (sy_mean=0.12, sy_sd=0.04),
    "SE_Coastal_Plain" => (sy_mean=0.10, sy_sd=0.04),
    "Floridan" => (sy_mean=0.001, sy_sd=0.0005),
    "NE_Glaciers" => (sy_mean=0.18, sy_sd=0.06),
    "Great_lakes" => (sy_mean=0.15, sy_sd=0.05),
    "N_Great_Plains" => (sy_mean=0.12, sy_sd=0.05),
    "Cambrian_Ordovician" => (sy_mean=0.0005, sy_sd=0.0003),
    "Piedmont_and_Blue_Ridge" => (sy_mean=0.01, sy_sd=0.005),
    "Ozark_Plateau" => (sy_mean=0.02, sy_sd=0.01),
    "Carrizo_Wilcox" => (sy_mean=0.0003, sy_sd=0.0002),
    "Trinity" => (sy_mean=0.0005, sy_sd=0.0003),
    "Edwards_Plateau" => (sy_mean=0.035, sy_sd=0.015),
    "Rio_Grande" => (sy_mean=0.15, sy_sd=0.05),
    "U_Colorado" => (sy_mean=0.08, sy_sd=0.04),
)

const DEFAULT_PRIOR = (sy_mean=0.10, sy_sd=0.05)

# =============================================================================
# Helper Functions
# =============================================================================

function sample_truncated_normal(μ, σ, n; lower=0.0001, upper=0.5)
    samples = Float64[]
    while length(samples) < n
        x = μ + σ * randn()
        if lower <= x <= upper
            push!(samples, x)
        end
    end
    return samples
end

safe_quantile(x, p) = isempty(filter(isfinite, x)) ? NaN : quantile(filter(isfinite, x), p)

function safe_mean(x)
    v = filter(v -> !ismissing(v) && isfinite(v), x)
    isempty(v) ? NaN : mean(v)
end

function safe_std(x)
    v = filter(v -> !ismissing(v) && isfinite(v), x)
    length(v) < 2 ? NaN : std(v)
end

function weighted_mean(values, weights)
    valid = [!ismissing(v) && isfinite(v) && !ismissing(w) && isfinite(w) && w > 0
             for (v, w) in zip(values, weights)]
    sum(valid) == 0 && return NaN
    v, w = values[valid], weights[valid]
    return sum(v .* w) / sum(w)
end

# =============================================================================
# Baseline Recomputation
# =============================================================================

function compute_baselines_for_period(df::DataFrame, baseline_start::Int, baseline_end::Int)
    """Compute well baselines for a specific period"""

    baseline_period = baseline_start:baseline_end
    baseline_data = df[in.(df.Year, Ref(Set(baseline_period))), :]

    if nrow(baseline_data) == 0
        return DataFrame()
    end

    # Aquifer-wide prior
    all_depths = filter(x -> !ismissing(x) && isfinite(x), baseline_data.DepthToWater_m)
    if isempty(all_depths)
        all_depths = filter(x -> !ismissing(x) && isfinite(x), df.DepthToWater_m)
    end
    aquifer_prior_mean = isempty(all_depths) ? 50.0 : median(all_depths)
    aquifer_prior_std = isempty(all_depths) ? 30.0 : mad(all_depths, normalize=true)

    # Bayesian baseline per well
    well_baselines = combine(groupby(baseline_data, :WellID)) do gdf
        depths = filter(x -> !ismissing(x) && isfinite(x), gdf.DepthToWater_m)
        years = gdf.Year[.!ismissing.(gdf.DepthToWater_m)]
        n_obs = length(depths)

        if n_obs < 2
            return (BaselineMean_m = missing, BaselineStd_m = missing, BaselineN = 0)
        end

        data_mean = mean(depths)
        data_std = n_obs > 1 ? std(depths) : aquifer_prior_std
        data_var = data_std^2 / n_obs

        age_factor = max(0.5, 1.0 - (baseline_end - mean(years)) / 50.0)
        n_factor = min(1.0, n_obs / 10.0)
        prior_var = (aquifer_prior_std^2) / (age_factor * n_factor)

        prior_precision = 1.0 / prior_var
        data_precision = 1.0 / data_var

        posterior_precision = prior_precision + data_precision
        posterior_mean = (prior_precision * aquifer_prior_mean + data_precision * data_mean) / posterior_precision
        posterior_std = sqrt(1.0 / posterior_precision)

        (BaselineMean_m = posterior_mean, BaselineStd_m = posterior_std, BaselineN = n_obs)
    end

    well_baselines = well_baselines[.!ismissing.(well_baselines.BaselineMean_m), :]
    return well_baselines
end

function compute_anomalies_and_storage(df::DataFrame, baselines::DataFrame,
                                       area_km2::Float64, sy_prior::NamedTuple)
    """Compute anomalies and storage for given baselines"""

    baseline_lookup = Dict(row.WellID => row.BaselineMean_m for row in eachrow(baselines))

    # Filter to wells with baselines
    df_anom = df[haskey.(Ref(baseline_lookup), df.WellID), :]
    if nrow(df_anom) == 0
        return (anomaly=NaN, storage_mean=NaN, storage_sd=NaN, n_wells=0)
    end

    df_anom = copy(df_anom)
    df_anom.Baseline_m = [baseline_lookup[wid] for wid in df_anom.WellID]
    df_anom.Anomaly_m = df_anom.Baseline_m .- df_anom.DepthToWater_m

    # Recent anomaly (2020+)
    recent = df_anom[df_anom.Year .>= 2020, :]
    if nrow(recent) == 0
        recent = df_anom[df_anom.Year .>= 2015, :]
    end

    if nrow(recent) == 0
        return (anomaly=NaN, storage_mean=NaN, storage_sd=NaN, n_wells=length(unique(df_anom.WellID)))
    end

    # IDW interpolation for spatial weighting
    valid = .!ismissing.(recent.Latitude) .& .!ismissing.(recent.Longitude) .&
            .!ismissing.(recent.Anomaly_m) .& isfinite.(recent.Anomaly_m)
    recent_valid = recent[valid, :]

    n_wells = length(unique(recent_valid.WellID))
    if n_wells == 0
        return (anomaly=NaN, storage_mean=NaN, storage_sd=NaN, n_wells=0)
    end

    z = Float64.(recent_valid.Anomaly_m)

    if n_wells >= 3
        x_m, y_m = project_to_aea(Float64.(recent_valid.Longitude), Float64.(recent_valid.Latitude))
        spacing_km = n_wells >= 5 ? choose_spacing_km(x_m, y_m) : 10.0
        grid_x, grid_y = make_grid_centers(x_m, y_m, spacing_km; buffer_km=spacing_km)

        if length(grid_x) > 0
            z_idw = idw_interpolate(x_m, y_m, z, grid_x, grid_y; idp=2, nmax=12)
            valid_idw = isfinite.(z_idw)
            idw_anomaly = sum(valid_idw) > 0 ? mean(z_idw[valid_idw]) : mean(z)
        else
            idw_anomaly = mean(z)
        end
    else
        idw_anomaly = mean(z)
    end

    anomaly_std = length(z) > 1 ? std(z) * 0.1 : abs(idw_anomaly) * 0.1

    # Monte Carlo storage
    sy_samples = sample_truncated_normal(sy_prior.sy_mean, sy_prior.sy_sd, N_SAMPLES)
    anomaly_samples = idw_anomaly .+ anomaly_std .* randn(N_SAMPLES)
    storage_samples = sy_samples .* area_km2 .* anomaly_samples ./ 1000.0

    return (
        anomaly = idw_anomaly,
        storage_mean = mean(storage_samples),
        storage_sd = std(storage_samples),
        n_wells = n_wells
    )
end

# =============================================================================
# Load Areas
# =============================================================================

function load_aquifer_areas()
    area_file = joinpath(INPUT_DIR, "aquifer_areas.csv")
    areas = CSV.read(area_file, DataFrame)
    return Dict(row.Aquifer_ID => row.Area_km2 for row in eachrow(areas))
end

# =============================================================================
# Main Analysis
# =============================================================================

function main()
    println("=" ^ 80)
    println("BASELINE PERIOD SENSITIVITY ANALYSIS")
    println("=" ^ 80)
    println()

    areas = load_aquifer_areas()
    println("Loaded $(length(areas)) aquifer areas")

    aquifer_dirs = filter(isdir, [joinpath(INPUT_DIR, d) for d in readdir(INPUT_DIR)])
    aquifer_dirs = filter(d -> isfile(joinpath(d, "anomalies.csv")), aquifer_dirs)
    println("Found $(length(aquifer_dirs)) aquifers with data")
    println()

    # Store all results
    all_results = DataFrame()
    aquifer_comparison = DataFrame()

    for (bp_start, bp_end, bp_name) in BASELINE_PERIODS
        println("\n" * "=" ^ 80)
        println("BASELINE PERIOD: $bp_name ($bp_start-$bp_end)")
        println("=" ^ 80)

        period_results = DataFrame()

        for (i, aq_dir) in enumerate(aquifer_dirs)
            aquifer = basename(aq_dir)
            area_km2 = get(areas, aquifer, NaN)

            isnan(area_km2) && continue

            anomaly_file = joinpath(aq_dir, "anomalies.csv")
            df = CSV.read(anomaly_file, DataFrame)

            nrow(df) < 20 && continue

            # Get Sy prior
            prior = get(SY_PRIORS, aquifer, DEFAULT_PRIOR)

            # Compute baselines for this period
            baselines = compute_baselines_for_period(df, bp_start, bp_end)

            if nrow(baselines) < 5
                continue
            end

            # Compute storage
            result = compute_anomalies_and_storage(df, baselines, area_km2, prior)

            if isnan(result.storage_mean)
                continue
            end

            push!(period_results, (
                Aquifer = aquifer,
                Baseline_Period = "$bp_start-$bp_end",
                Area_km2 = round(area_km2, digits=0),
                N_Wells = result.n_wells,
                Anomaly_m = round(result.anomaly, digits=2),
                Storage_km3 = round(result.storage_mean, digits=1),
                Storage_SD_km3 = round(result.storage_sd, digits=1)
            ))
        end

        # Summary for this period
        if nrow(period_results) > 0
            total_storage = sum(period_results.Storage_km3)
            total_sd = sqrt(sum(period_results.Storage_SD_km3 .^ 2))

            println("\n  Aquifers: $(nrow(period_results))")
            println("  Total storage: $(round(total_storage, digits=0)) ± $(round(total_sd, digits=0)) km³")

            push!(aquifer_comparison, (
                Baseline = "$bp_start-$bp_end",
                Label = bp_name,
                N_Aquifers = nrow(period_results),
                Total_Storage_km3 = round(total_storage, digits=0),
                Total_SD_km3 = round(total_sd, digits=0),
                CI_Lower_km3 = round(total_storage - 1.96*total_sd, digits=0),
                CI_Upper_km3 = round(total_storage + 1.96*total_sd, digits=0)
            ))

            append!(all_results, period_results)
        end
    end

    # =============================================================================
    # Save Results
    # =============================================================================

    println("\n" * "=" ^ 80)
    println("SENSITIVITY SUMMARY")
    println("=" ^ 80)

    println("\n" * "-" ^ 60)
    println("Total CONUS Storage by Baseline Period")
    println("-" ^ 60)

    for row in eachrow(aquifer_comparison)
        println("  $(rpad(row.Label, 25)) $(lpad(Int(row.Total_Storage_km3), 7)) ± $(lpad(Int(row.Total_SD_km3), 4)) km³")
    end

    # Compute sensitivity range
    storages = aquifer_comparison.Total_Storage_km3
    range_km3 = maximum(storages) - minimum(storages)
    mean_storage = mean(storages)
    range_pct = 100 * range_km3 / abs(mean_storage)

    println("\n" * "-" ^ 60)
    println("Sensitivity Metrics")
    println("-" ^ 60)
    println("  Storage range:      $(round(range_km3, digits=0)) km³")
    println("  Range as % of mean: $(round(range_pct, digits=1))%")
    println("  Mean across periods: $(round(mean_storage, digits=0)) km³")
    println("  Std across periods:  $(round(std(storages), digits=0)) km³")

    # Save files
    CSV.write(joinpath(OUTPUT_DIR, "baseline_sensitivity_all.csv"), all_results)
    CSV.write(joinpath(OUTPUT_DIR, "baseline_sensitivity_summary.csv"), aquifer_comparison)

    # Create per-aquifer comparison
    println("\n" * "-" ^ 60)
    println("Per-Aquifer Sensitivity (Top 10 most sensitive)")
    println("-" ^ 60)

    # Pivot to wide format for comparison
    aquifer_wide = DataFrame()
    for aquifer in unique(all_results.Aquifer)
        aq_data = all_results[all_results.Aquifer .== aquifer, :]
        if nrow(aq_data) == length(BASELINE_PERIODS)
            storages = aq_data.Storage_km3
            row = (
                Aquifer = aquifer,
                Storage_1930_1945 = storages[aq_data.Baseline_Period .== "1930-1945"][1],
                Storage_1935_1950 = storages[aq_data.Baseline_Period .== "1935-1950"][1],
                Storage_1940_1955 = storages[aq_data.Baseline_Period .== "1940-1955"][1],
                Storage_1945_1960 = storages[aq_data.Baseline_Period .== "1945-1960"][1],
                Storage_1950_1965 = storages[aq_data.Baseline_Period .== "1950-1965"][1],
                Range_km3 = maximum(storages) - minimum(storages),
                Range_Pct = 100 * (maximum(storages) - minimum(storages)) / abs(mean(storages))
            )
            push!(aquifer_wide, row)
        end
    end

    if nrow(aquifer_wide) > 0
        sort!(aquifer_wide, :Range_km3, rev=true)

        for row in eachrow(first(aquifer_wide, 10))
            println("  $(rpad(row.Aquifer, 35)) Range: $(lpad(round(Int, row.Range_km3), 5)) km³ ($(round(row.Range_Pct, digits=1))%)")
        end

        CSV.write(joinpath(OUTPUT_DIR, "baseline_sensitivity_aquifer_comparison.csv"), aquifer_wide)
    end

    println("\n" * "=" ^ 80)
    println("OUTPUT FILES:")
    println("  $(joinpath(OUTPUT_DIR, "baseline_sensitivity_summary.csv"))")
    println("  $(joinpath(OUTPUT_DIR, "baseline_sensitivity_all.csv"))")
    println("  $(joinpath(OUTPUT_DIR, "baseline_sensitivity_aquifer_comparison.csv"))")
    println("=" ^ 80)

    return aquifer_comparison
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

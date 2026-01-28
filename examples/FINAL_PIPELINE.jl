#!/usr/bin/env julia
# =============================================================================
# FINAL GROUNDWATER STORAGE ANALYSIS PIPELINE
# =============================================================================
#
# Complete A-Z analysis for CONUS groundwater storage estimation
#
# Author: Dr. Ashraf Rateb, Bureau of Economic Geology, UT Austin
# Date: January 2026
#
# =============================================================================

using Pkg
Pkg.activate("/data2/GWProject")

using DataFrames, CSV, Statistics, Dates, Random
using StatsBase: quantile, std, mad, iqr

Random.seed!(42)

# Load spatial methods
include("/data2/GWProject/src/spatial_methods.jl")
using .SpatialMethods

# =============================================================================
# CONFIGURATION
# =============================================================================

const CONFIG = (
    # Paths
    input_csv = "/data2/GWProject/data/CONUS/GW_CONUS_Harmonized_20251208.csv",
    aquifer_mapping = "/data2/GWProject/outputs/Phase1/well_aquifer_mapping_DEDUPED.csv",
    output_base = "/data2/GWProject/outputs",

    # Temporal
    baseline_start = 1940,
    baseline_end = 1955,
    year_min = 1940,
    year_max = 2025,
    valid_seasons = ["Winter", "Spring"],

    # Quality
    min_obs_per_well = 5,
    min_years_baseline = 2,
    min_wells_aquifer = 5,
    outlier_zscore = 3.5,

    # Spatial
    idw_power = 2,
    idw_neighbors = 12,
    min_wells_idw = 3,

    # Monte Carlo
    n_samples = 10000
)

println("=" ^ 80)
println("GROUNDWATER STORAGE ANALYSIS - FINAL PIPELINE")
println("=" ^ 80)
println("Baseline period: $(CONFIG.baseline_start)-$(CONFIG.baseline_end)")
println("Monte Carlo samples: $(CONFIG.n_samples)")
println("Output: $(CONFIG.output_base)")
println()

# =============================================================================
# SPECIFIC YIELD PRIORS (Literature-based)
# =============================================================================

const SY_PRIORS = Dict(
    "N_High_Plains" => (sy=0.15, sd=0.04, conf=:unconfined, ref="Gutentag 1984"),
    "C_S_High_Plains" => (sy=0.15, sd=0.04, conf=:unconfined, ref="Gutentag 1984"),
    "High_Plains" => (sy=0.15, sd=0.04, conf=:unconfined, ref="Gutentag 1984"),
    "Central_Valley" => (sy=0.12, sd=0.04, conf=:mixed, ref="Faunt 2009"),
    "Basin_and_Range_Basin_Fill" => (sy=0.18, sd=0.04, conf=:mixed, ref="Pool & Coes 1999"),
    "Arizona_Alluvials" => (sy=0.18, sd=0.05, conf=:unconfined, ref="Pool & Coes 1999"),
    "Mississippi_River_Valley_Alluvial" => (sy=0.22, sd=0.05, conf=:unconfined, ref="USGS MERAS"),
    "Mississippi_Embayment" => (sy=0.12, sd=0.04, conf=:mixed, ref="Clark & Hart 2009"),
    "Snake" => (sy=0.08, sd=0.04, conf=:unconfined, ref="Whitehead 1992"),
    "Columbia_Plateau_Basin_fill" => (sy=0.12, sd=0.05, conf=:unconfined, ref="USGS Columbia"),
    "Columbia_Plateau_Basaltic_rock" => (sy=0.05, sd=0.03, conf=:mixed, ref="USGS Columbia"),
    "Coastal_Lowland" => (sy=0.08, sd=0.04, conf=:mixed, ref="USGS Gulf Coast"),
    "Texas_Gulf_Coast" => (sy=0.06, sd=0.03, conf=:confined, ref="USGS Gulf Coast"),
    "N_Atlantic_Coastal_Plain" => (sy=0.12, sd=0.04, conf=:mixed, ref="USGS Atlantic"),
    "SE_Coastal_Plain" => (sy=0.10, sd=0.04, conf=:mixed, ref="USGS Coastal Plain"),
    "Floridan" => (sy=0.001, sd=0.0005, conf=:confined, ref="Miller 1986"),
    "NE_Glaciers" => (sy=0.18, sd=0.06, conf=:unconfined, ref="USGS glacial"),
    "Great_lakes" => (sy=0.15, sd=0.05, conf=:mixed, ref="USGS Great Lakes"),
    "N_Great_Plains" => (sy=0.12, sd=0.05, conf=:unconfined, ref="USGS N Great Plains"),
    "Cambrian_Ordovician" => (sy=0.0005, sd=0.0003, conf=:confined, ref="Young 1992"),
    "Piedmont_and_Blue_Ridge" => (sy=0.01, sd=0.005, conf=:mixed, ref="Daniel 1989"),
    "Ozark_Plateau" => (sy=0.02, sd=0.01, conf=:confined, ref="Imes & Emmett 1994"),
    "Carrizo_Wilcox" => (sy=0.0003, sd=0.0002, conf=:confined, ref="Mace 2000"),
    "Trinity" => (sy=0.0005, sd=0.0003, conf=:confined, ref="TWDB GAM"),
    "Edwards_Plateau" => (sy=0.035, sd=0.015, conf=:mixed, ref="Maclay & Small 1986"),
    "Rio_Grande" => (sy=0.15, sd=0.05, conf=:unconfined, ref="USGS Rio Grande"),
    "U_Colorado" => (sy=0.08, sd=0.04, conf=:mixed, ref="USGS Upper Colorado"),
)

const DEFAULT_SY = (sy=0.10, sd=0.05, conf=:mixed, ref="Default")

# =============================================================================
# AQUIFER AREAS (km², EPSG:5070)
# =============================================================================

const AQUIFER_AREAS = Dict(
    "N_High_Plains" => 248320.0,
    "C_S_High_Plains" => 203400.0,
    "High_Plains" => 505486.0,
    "Central_Valley" => 62421.0,
    "Basin_and_Range_Basin_Fill" => 429865.0,
    "Arizona_Alluvials" => 209945.0,
    "Mississippi_River_Valley_Alluvial" => 94635.0,
    "Mississippi_Embayment" => 199098.0,
    "Snake" => 186489.0,
    "Columbia_Plateau_Basin_fill" => 62683.0,
    "Columbia_Plateau_Basaltic_rock" => 110957.0,
    "Coastal_Lowland" => 143031.0,
    "Texas_Gulf_Coast" => 126522.0,
    "N_Atlantic_Coastal_Plain" => 134847.0,
    "SE_Coastal_Plain" => 125424.0,
    "Floridan" => 214839.0,
    "NE_Glaciers" => 244101.0,
    "Great_lakes" => 653488.0,
    "N_Great_Plains" => 670433.0,
    "Cambrian_Ordovician" => 152207.0,
    "Piedmont_and_Blue_Ridge" => 376474.0,
    "Ozark_Plateau" => 134633.0,
    "Carrizo_Wilcox" => 104887.0,
    "Trinity" => 98266.0,
    "Edwards_Plateau" => 98126.0,
    "Rio_Grande" => 112813.0,
    "U_Colorado" => 369858.0,
)

# =============================================================================
# HELPER FUNCTIONS
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

safe_mean(x) = begin
    v = filter(v -> !ismissing(v) && isfinite(v), x)
    isempty(v) ? NaN : mean(v)
end

safe_std(x) = begin
    v = filter(v -> !ismissing(v) && isfinite(v), x)
    length(v) < 2 ? NaN : std(v)
end

safe_quantile(x, p) = begin
    v = filter(isfinite, x)
    isempty(v) ? NaN : quantile(v, p)
end

# =============================================================================
# PHASE 1: DATA LOADING AND CLEANING
# =============================================================================

println("=" ^ 80)
println("PHASE 1: DATA LOADING AND CLEANING")
println("=" ^ 80)

# Load harmonized data
println("Loading harmonized data...")
df_raw = CSV.read(CONFIG.input_csv, DataFrame)
df_raw.WellID = string.(df_raw.WellID)
println("  Raw records: $(nrow(df_raw))")

# Filter years
df = df_raw[df_raw.Year .>= CONFIG.year_min, :]
df = df[df.Year .<= CONFIG.year_max, :]
println("  After year filter ($(CONFIG.year_min)-$(CONFIG.year_max)): $(nrow(df))")

# Filter CONUS bounds
df = df[(df.Latitude .>= 24.0) .& (df.Latitude .<= 50.0), :]
df = df[(df.Longitude .>= -125.0) .& (df.Longitude .<= -66.0), :]
println("  After CONUS bounds: $(nrow(df))")

# Filter valid depths
df = df[df.DepthToWater_m .> 0, :]
df = df[df.DepthToWater_m .< 500, :]  # Reasonable max depth
println("  After depth filter (0-500m): $(nrow(df))")

# Filter seasons (Winter/Spring = pre-irrigation)
if :Season in propertynames(df)
    df = df[in.(df.Season, Ref(Set(CONFIG.valid_seasons))), :]
    println("  After season filter (Winter/Spring): $(nrow(df))")
end

# Load aquifer mapping
println("\nLoading aquifer mapping...")

# Function to standardize aquifer names
function standardize_aquifer_name(name::AbstractString)
    s = replace(name, "+" => "_", "-" => "_", " " => "_")
    s = replace(s, r"_+" => "_")  # collapse multiple underscores
    return s
end

if !isfile(CONFIG.aquifer_mapping)
    println("  Using existing Aquifer column from harmonized data")
    df.Aquifer_Clean = standardize_aquifer_name.(coalesce.(df.Aquifer, "Unknown"))
else
    mapping = CSV.read(CONFIG.aquifer_mapping, DataFrame)
    mapping.WellID = string.(mapping.WellID)
    # Standardize aquifer names in mapping
    mapping.Aquifer_Std = standardize_aquifer_name.(mapping.Aquifer_EPA)
    mapping_dict = Dict(row.WellID => row.Aquifer_Std for row in eachrow(mapping))
    df.Aquifer_Clean = [get(mapping_dict, wid, "Unknown") for wid in df.WellID]
    println("  Mapped $(sum(df.Aquifer_Clean .!= "Unknown")) wells to aquifers")
end

# Remove unknown aquifers
df = df[df.Aquifer_Clean .!= "Unknown", :]
println("  After removing Unknown: $(nrow(df))")

# Get unique aquifers
aquifers = sort(unique(df.Aquifer_Clean))
println("  Unique aquifers in data: $(length(aquifers))")
aquifers = filter(a -> haskey(AQUIFER_AREAS, a), aquifers)
println("  Aquifers with area data: $(length(aquifers))")

# Save Phase 1
phase1_dir = joinpath(CONFIG.output_base, "Phase1_Clean")
mkpath(phase1_dir)
CSV.write(joinpath(phase1_dir, "GW_Cleaned.csv"), df)
println("Saved: $(joinpath(phase1_dir, "GW_Cleaned.csv"))")

# =============================================================================
# PHASE 2: BASELINES AND ANOMALIES
# =============================================================================

println("\n" * "=" ^ 80)
println("PHASE 2: BASELINES AND ANOMALIES")
println("=" ^ 80)

phase2_dir = joinpath(CONFIG.output_base, "Phase2_Baselines")
mkpath(phase2_dir)

all_baselines = DataFrame()
all_anomalies = DataFrame()

for aquifer in aquifers
    df_aq = df[df.Aquifer_Clean .== aquifer, :]
    n_wells = length(unique(df_aq.WellID))

    if n_wells < CONFIG.min_wells_aquifer
        continue
    end

    # Baseline period data
    baseline_data = df_aq[(df_aq.Year .>= CONFIG.baseline_start) .&
                          (df_aq.Year .<= CONFIG.baseline_end), :]

    if nrow(baseline_data) == 0
        continue
    end

    # Compute per-well baselines
    well_baselines = combine(groupby(baseline_data, :WellID)) do gdf
        depths = filter(x -> !ismissing(x) && isfinite(x), gdf.DepthToWater_m)
        n_years = length(unique(gdf.Year))

        if n_years < CONFIG.min_years_baseline
            return (Baseline_m=missing, Baseline_SD=missing, N_Years=0)
        end

        (Baseline_m=mean(depths), Baseline_SD=std(depths), N_Years=n_years)
    end

    well_baselines = well_baselines[.!ismissing.(well_baselines.Baseline_m), :]
    well_baselines.Aquifer .= aquifer

    if nrow(well_baselines) < CONFIG.min_wells_aquifer
        continue
    end

    # Compute anomalies
    baseline_lookup = Dict(row.WellID => row.Baseline_m for row in eachrow(well_baselines))

    df_anom = df_aq[haskey.(Ref(baseline_lookup), df_aq.WellID), :]
    df_anom = copy(df_anom)
    df_anom.Baseline_m = [baseline_lookup[wid] for wid in df_anom.WellID]
    # GRACE convention: Anomaly = Baseline - Current (negative = depletion)
    df_anom.Anomaly_m = df_anom.Baseline_m .- df_anom.DepthToWater_m

    append!(all_baselines, well_baselines)
    append!(all_anomalies, df_anom)

    print("  $(rpad(aquifer, 35)) $(lpad(nrow(well_baselines), 5)) wells\n")
end

CSV.write(joinpath(phase2_dir, "well_baselines.csv"), all_baselines)
CSV.write(joinpath(phase2_dir, "anomalies.csv"), all_anomalies)

println("\nTotal wells with baseline: $(nrow(all_baselines))")
println("Total records with anomalies: $(nrow(all_anomalies))")

# =============================================================================
# PHASE 3: STORAGE ESTIMATION
# =============================================================================

println("\n" * "=" ^ 80)
println("PHASE 3: BAYESIAN STORAGE ESTIMATION")
println("=" ^ 80)

phase3_dir = joinpath(CONFIG.output_base, "Phase3_Storage")
mkpath(phase3_dir)
mkpath(joinpath(phase3_dir, "annual"))

storage_results = DataFrame()

for aquifer in unique(all_anomalies.Aquifer_Clean)
    df_aq = all_anomalies[all_anomalies.Aquifer_Clean .== aquifer, :]
    area_km2 = get(AQUIFER_AREAS, aquifer, NaN)
    sy_prior = get(SY_PRIORS, aquifer, DEFAULT_SY)

    if isnan(area_km2) || nrow(df_aq) < 100
        continue
    end

    # Compute annual anomalies with IDW
    years = sort(unique(df_aq.Year))
    annual = DataFrame()

    for year in years
        df_year = df_aq[df_aq.Year .== year, :]
        valid = .!ismissing.(df_year.Anomaly_m) .& isfinite.(df_year.Anomaly_m)
        df_valid = df_year[valid, :]

        n_wells = length(unique(df_valid.WellID))
        n_wells == 0 && continue

        z = Float64.(df_valid.Anomaly_m)
        simple_mean = mean(z)

        # IDW interpolation
        if n_wells >= CONFIG.min_wells_idw
            x_m, y_m = project_to_aea(Float64.(df_valid.Longitude), Float64.(df_valid.Latitude))
            spacing_km = n_wells >= 5 ? choose_spacing_km(x_m, y_m) : 10.0
            grid_x, grid_y = make_grid_centers(x_m, y_m, spacing_km; buffer_km=spacing_km)

            if length(grid_x) > 0
                z_idw = idw_interpolate(x_m, y_m, z, grid_x, grid_y;
                                        idp=CONFIG.idw_power, nmax=CONFIG.idw_neighbors)
                valid_idw = isfinite.(z_idw)
                idw_mean = sum(valid_idw) > 0 ? mean(z_idw[valid_idw]) : simple_mean
                method = "IDW"
            else
                idw_mean = simple_mean
                method = "Simple"
            end
        else
            idw_mean = simple_mean
            method = "Simple"
        end

        push!(annual, (Year=year, Simple_Mean=simple_mean, IDW_Mean=idw_mean,
                       N_Wells=n_wells, Method=method))
    end

    if nrow(annual) < 10
        continue
    end

    # Get recent anomaly (2020+ or 2015+)
    recent = annual[annual.Year .>= 2020, :]
    if nrow(recent) == 0
        recent = annual[annual.Year .>= 2015, :]
    end
    nrow(recent) == 0 && continue

    idw_anomaly = safe_mean(recent.IDW_Mean)
    simple_anomaly = safe_mean(recent.Simple_Mean)
    anomaly_std = nrow(recent) > 1 ? safe_std(recent.IDW_Mean) : abs(idw_anomaly) * 0.1

    isnan(idw_anomaly) && continue

    # Linear trend
    valid_annual = annual[isfinite.(annual.IDW_Mean), :]
    if nrow(valid_annual) >= 10
        yrs = Float64.(valid_annual.Year)
        vals = Float64.(valid_annual.IDW_Mean)
        x_mean, y_mean = mean(yrs), mean(vals)
        Sxx = sum((yrs .- x_mean).^2)
        Sxy = sum((yrs .- x_mean) .* (vals .- y_mean))
        slope = Sxy / Sxx
        y_pred = y_mean .+ slope .* (yrs .- x_mean)
        ss_res = sum((vals .- y_pred).^2)
        ss_tot = sum((vals .- y_mean).^2)
        r2 = ss_tot > 0 ? 1 - ss_res / ss_tot : NaN
        trend_m_dec = slope * 10
    else
        trend_m_dec, r2 = NaN, NaN
    end

    # Monte Carlo storage estimation
    sy_samples = sample_truncated_normal(sy_prior.sy, sy_prior.sd, CONFIG.n_samples)
    anomaly_samples = idw_anomaly .+ anomaly_std .* randn(CONFIG.n_samples)
    storage_samples = sy_samples .* area_km2 .* anomaly_samples ./ 1000.0

    storage_mean = mean(storage_samples)
    storage_sd = std(storage_samples)

    # Rate samples
    if !isnan(trend_m_dec)
        trend_se = abs(trend_m_dec) * 0.1
        trend_samples = trend_m_dec .+ trend_se .* randn(CONFIG.n_samples)
        rate_samples = sy_samples .* area_km2 .* trend_samples ./ 1000.0
        rate_mean = mean(rate_samples)
        rate_sd = std(rate_samples)
    else
        rate_mean, rate_sd = NaN, NaN
    end

    spatial_method = "IDW" in recent.Method ? "IDW" : "Simple"

    push!(storage_results, (
        Aquifer = aquifer,
        Area_km2 = area_km2,
        Confinement = String(sy_prior.conf),
        Sy_Mean = sy_prior.sy,
        Sy_SD = sy_prior.sd,
        Sy_Reference = sy_prior.ref,
        N_Years = nrow(annual),
        Year_Start = minimum(annual.Year),
        Year_End = maximum(annual.Year),
        N_Wells = round(Int, mean(annual.N_Wells)),
        Spatial_Method = spatial_method,
        Simple_Anomaly_m = round(simple_anomaly, digits=2),
        IDW_Anomaly_m = round(idw_anomaly, digits=2),
        Spatial_Correction_m = round(idw_anomaly - simple_anomaly, digits=2),
        Trend_m_decade = round(trend_m_dec, digits=3),
        Trend_R2 = round(r2, digits=3),
        Storage_km3 = round(storage_mean, digits=1),
        Storage_SD_km3 = round(storage_sd, digits=1),
        Storage_Q05_km3 = round(safe_quantile(storage_samples, 0.05), digits=1),
        Storage_Q95_km3 = round(safe_quantile(storage_samples, 0.95), digits=1),
        Rate_km3_decade = round(rate_mean, digits=1),
        Rate_SD_km3_decade = round(rate_sd, digits=1)
    ))

    # Save annual time series
    CSV.write(joinpath(phase3_dir, "annual", "$(aquifer)_annual.csv"), annual)

    symbol = storage_mean < 0 ? "▼" : "▲"
    println("  $symbol $(rpad(aquifer, 30)) $(lpad(round(Int, storage_mean), 6)) ± $(lpad(round(Int, storage_sd), 4)) km³")
end

# Sort by storage (most depleted first)
sort!(storage_results, :Storage_km3)

CSV.write(joinpath(phase3_dir, "STORAGE_RESULTS.csv"), storage_results)

# Total CONUS
total_storage = sum(storage_results.Storage_km3)
total_sd = sqrt(sum(storage_results.Storage_SD_km3 .^ 2))

println("\n" * "-" ^ 60)
println("TOTAL CONUS: $(round(Int, total_storage)) ± $(round(Int, total_sd)) km³")
println("-" ^ 60)

# =============================================================================
# PHASE 4: BASELINE SENSITIVITY
# =============================================================================

println("\n" * "=" ^ 80)
println("PHASE 4: BASELINE SENSITIVITY ANALYSIS")
println("=" ^ 80)

phase4_dir = joinpath(CONFIG.output_base, "Phase4_Sensitivity")
mkpath(phase4_dir)

baseline_periods = [
    (1930, 1945, "1930-1945"),
    (1935, 1950, "1935-1950"),
    (1940, 1955, "1940-1955"),
    (1945, 1960, "1945-1960"),
    (1950, 1965, "1950-1965")
]

sensitivity_summary = DataFrame()

for (bp_start, bp_end, bp_label) in baseline_periods
    # Recompute baselines for this period
    period_storage = 0.0
    period_sd_sq = 0.0
    n_aquifers = 0

    for aquifer in unique(df.Aquifer_Clean)
        df_aq = df[df.Aquifer_Clean .== aquifer, :]
        area_km2 = get(AQUIFER_AREAS, aquifer, NaN)
        sy_prior = get(SY_PRIORS, aquifer, DEFAULT_SY)

        isnan(area_km2) && continue

        # Baseline for this period
        baseline_data = df_aq[(df_aq.Year .>= bp_start) .& (df_aq.Year .<= bp_end), :]
        nrow(baseline_data) < 10 && continue

        well_baselines = combine(groupby(baseline_data, :WellID)) do gdf
            depths = filter(x -> !ismissing(x) && isfinite(x), gdf.DepthToWater_m)
            length(depths) < 2 && return (Baseline_m=missing,)
            (Baseline_m=mean(depths),)
        end
        well_baselines = well_baselines[.!ismissing.(well_baselines.Baseline_m), :]
        nrow(well_baselines) < 5 && continue

        # Recent data
        recent_data = df_aq[df_aq.Year .>= 2015, :]
        nrow(recent_data) < 5 && continue

        baseline_lookup = Dict(row.WellID => row.Baseline_m for row in eachrow(well_baselines))
        recent_with_bl = recent_data[haskey.(Ref(baseline_lookup), recent_data.WellID), :]
        nrow(recent_with_bl) < 5 && continue

        anomalies = [baseline_lookup[wid] - d for (wid, d) in
                     zip(recent_with_bl.WellID, recent_with_bl.DepthToWater_m)]
        anomaly_mean = mean(filter(isfinite, anomalies))
        isnan(anomaly_mean) && continue

        # Storage
        storage = sy_prior.sy * area_km2 * anomaly_mean / 1000.0
        storage_sd = sy_prior.sd * area_km2 * abs(anomaly_mean) / 1000.0

        period_storage += storage
        period_sd_sq += storage_sd^2
        n_aquifers += 1
    end

    period_sd = sqrt(period_sd_sq)

    push!(sensitivity_summary, (
        Baseline_Period = bp_label,
        N_Aquifers = n_aquifers,
        Total_Storage_km3 = round(period_storage, digits=0),
        Total_SD_km3 = round(period_sd, digits=0)
    ))

    println("  $(rpad(bp_label, 15)) $(lpad(round(Int, period_storage), 7)) ± $(lpad(round(Int, period_sd), 4)) km³  ($n_aquifers aquifers)")
end

CSV.write(joinpath(phase4_dir, "baseline_sensitivity.csv"), sensitivity_summary)

# Sensitivity range
range_km3 = maximum(sensitivity_summary.Total_Storage_km3) - minimum(sensitivity_summary.Total_Storage_km3)
mean_storage = mean(sensitivity_summary.Total_Storage_km3)
println("\n  Sensitivity range: $(round(Int, range_km3)) km³ ($(round(100*range_km3/abs(mean_storage), digits=1))% of mean)")

# =============================================================================
# PHASE 5: TREND BREAKPOINTS
# =============================================================================

println("\n" * "=" ^ 80)
println("PHASE 5: TREND BREAKPOINT ANALYSIS")
println("=" ^ 80)

phase5_dir = joinpath(CONFIG.output_base, "Phase5_Trends")
mkpath(phase5_dir)

breakpoint_results = DataFrame()

for row in eachrow(storage_results)
    aquifer = row.Aquifer
    annual_file = joinpath(phase3_dir, "annual", "$(aquifer)_annual.csv")

    !isfile(annual_file) && continue

    annual = CSV.read(annual_file, DataFrame)
    valid = isfinite.(annual.IDW_Mean)
    sum(valid) < 20 && continue

    years = Float64.(annual.Year[valid])
    values = Float64.(annual.IDW_Mean[valid])

    # Single linear trend
    x_mean, y_mean = mean(years), mean(values)
    Sxx = sum((years .- x_mean).^2)
    Sxy = sum((years .- x_mean) .* (values .- y_mean))
    slope_single = Sxy / Sxx
    pred_single = y_mean .+ slope_single .* (years .- x_mean)
    rss_single = sum((values .- pred_single).^2)
    tss = sum((values .- y_mean).^2)
    r2_single = 1 - rss_single / tss

    # Find best breakpoint
    best_bp, best_rss = 0, Inf
    best_slope1, best_slope2 = NaN, NaN

    for bp in (Int(minimum(years))+10):5:(Int(maximum(years))-10)
        m1 = years .<= bp
        m2 = years .> bp

        (sum(m1) < 5 || sum(m2) < 5) && continue

        # Fit segments
        y1, v1 = years[m1], values[m1]
        y2, v2 = years[m2], values[m2]

        slope1 = sum((y1 .- mean(y1)) .* (v1 .- mean(v1))) / sum((y1 .- mean(y1)).^2)
        int1 = mean(v1) - slope1 * mean(y1)

        slope2 = sum((y2 .- mean(y2)) .* (v2 .- mean(v2))) / sum((y2 .- mean(y2)).^2)
        int2 = mean(v2) - slope2 * mean(y2)

        rss = sum((v1 .- (int1 .+ slope1 .* y1)).^2) + sum((v2 .- (int2 .+ slope2 .* y2)).^2)

        if rss < best_rss
            best_rss = rss
            best_bp = bp
            best_slope1 = slope1 * 10
            best_slope2 = slope2 * 10
        end
    end

    r2_piecewise = 1 - best_rss / tss
    improvement = r2_piecewise - r2_single
    acceleration = best_slope2 - best_slope1

    # Classify pattern
    pattern = if improvement < 0.02
        "LINEAR"
    elseif acceleration < -0.5
        "ACCELERATING"
    elseif acceleration > 0.5
        "DECELERATING"
    elseif best_slope1 < 0 && best_slope2 > 0
        "RECOVERY"
    else
        "STABLE"
    end

    push!(breakpoint_results, (
        Aquifer = aquifer,
        Breakpoint = best_bp,
        Slope_Before = round(best_slope1, digits=2),
        Slope_After = round(best_slope2, digits=2),
        Acceleration = round(acceleration, digits=2),
        R2_Improvement = round(improvement, digits=3),
        Pattern = pattern
    ))

    if improvement >= 0.05
        println("  $(rpad(aquifer, 30)) $(best_bp): $(round(best_slope1, digits=1)) → $(round(best_slope2, digits=1)) m/dec  $pattern")
    end
end

CSV.write(joinpath(phase5_dir, "trend_breakpoints.csv"), breakpoint_results)

# Pattern summary
pattern_counts = combine(groupby(breakpoint_results, :Pattern), nrow => :Count)
println("\nTrend patterns:")
for row in eachrow(sort(pattern_counts, :Count, rev=true))
    println("  $(rpad(row.Pattern, 15)) $(row.Count) aquifers")
end

# =============================================================================
# FINAL SUMMARY
# =============================================================================

println("\n" * "=" ^ 80)
println("FINAL SUMMARY")
println("=" ^ 80)

final_dir = joinpath(CONFIG.output_base, "Final")
mkpath(final_dir)

# Copy main results to Final
cp(joinpath(phase3_dir, "STORAGE_RESULTS.csv"),
   joinpath(final_dir, "STORAGE_RESULTS.csv"), force=true)

# Create summary report
summary_text = """
# GROUNDWATER STORAGE ANALYSIS - FINAL RESULTS

**Date**: $(Dates.today())
**Author**: Dr. Ashraf Rateb, Bureau of Economic Geology, UT Austin

---

## Key Results

**Total CONUS Storage Change**: $(round(Int, total_storage)) ± $(round(Int, total_sd)) km³

**Period**: $(CONFIG.year_min)-$(CONFIG.year_max) relative to $(CONFIG.baseline_start)-$(CONFIG.baseline_end) baseline

**Aquifers Analyzed**: $(nrow(storage_results))

---

## Top 10 Storage Changes

| Rank | Aquifer | Storage (km³) | Rate (km³/dec) |
|------|---------|---------------|----------------|
"""

for (i, row) in enumerate(eachrow(first(storage_results, 10)))
    global summary_text *= "| $i | $(row.Aquifer) | $(row.Storage_km3) ± $(row.Storage_SD_km3) | $(row.Rate_km3_decade) |\n"
end

global summary_text *= """

---

## Baseline Sensitivity

| Period | Storage (km³) |
|--------|---------------|
"""

for row in eachrow(sensitivity_summary)
    global summary_text *= "| $(row.Baseline_Period) | $(row.Total_Storage_km3) ± $(row.Total_SD_km3) |\n"
end

summary_text *= """

**Sensitivity Range**: $(round(Int, range_km3)) km³ ($(round(100*range_km3/abs(mean_storage), digits=1))% of mean)

---

## Trend Patterns

"""

for row in eachrow(sort(pattern_counts, :Count, rev=true))
    global summary_text *= "- **$(row.Pattern)**: $(row.Count) aquifers\n"
end

summary_text *= """

---

## Output Files

```
$(CONFIG.output_base)/
├── Phase1_Clean/
│   └── GW_Cleaned.csv
├── Phase2_Baselines/
│   ├── well_baselines.csv
│   └── anomalies.csv
├── Phase3_Storage/
│   ├── STORAGE_RESULTS.csv
│   └── annual/*.csv
├── Phase4_Sensitivity/
│   └── baseline_sensitivity.csv
├── Phase5_Trends/
│   └── trend_breakpoints.csv
└── Final/
    ├── STORAGE_RESULTS.csv
    └── SUMMARY.md
```

---

*Generated by FINAL_PIPELINE.jl*
"""

open(joinpath(final_dir, "SUMMARY.md"), "w") do f
    write(f, summary_text)
end

println("\n" * "=" ^ 80)
println("PIPELINE COMPLETE")
println("=" ^ 80)
println("\nTotal CONUS: $(round(Int, total_storage)) ± $(round(Int, total_sd)) km³")
println("\nOutputs: $(CONFIG.output_base)/")
println("  - Phase1_Clean/GW_Cleaned.csv")
println("  - Phase2_Baselines/{well_baselines,anomalies}.csv")
println("  - Phase3_Storage/STORAGE_RESULTS.csv")
println("  - Phase4_Sensitivity/baseline_sensitivity.csv")
println("  - Phase5_Trends/trend_breakpoints.csv")
println("  - Final/SUMMARY.md")
println("=" ^ 80)

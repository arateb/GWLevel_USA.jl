# =============================================================================
# High-Level Analysis Functions
# =============================================================================

"""
    run_aquifer_analysis(df::DataFrame, aquifer::String, area_km2::Real;
                         cfg::GWConfig=DEFAULT_CONFIG) -> AquiferResult

Run complete analysis for a single aquifer.

# Arguments
- `df`: DataFrame with groundwater data (already filtered to aquifer)
- `aquifer`: Aquifer name
- `area_km2`: Aquifer area in km²
- `cfg`: Configuration parameters

# Returns
AquiferResult with storage estimates, trends, and annual data
"""
function run_aquifer_analysis(
    df::DataFrame,
    aquifer::String,
    area_km2::Real;
    cfg::GWConfig=DEFAULT_CONFIG
)
    # Get Sy prior
    sy_prior = get_sy_prior(aquifer)

    # Filter to valid data
    df = df[isfinite.(df.DepthToWater_m), :]

    if nrow(df) < cfg.min_records_per_aquifer
        @warn "Insufficient data for $aquifer: $(nrow(df)) records"
        return nothing
    end

    # Compute baselines
    baselines = compute_baseline_all(df; cfg=cfg)

    if nrow(baselines) < cfg.min_wells_per_aquifer
        @warn "Insufficient wells with baseline for $aquifer: $(nrow(baselines))"
        return nothing
    end

    # Add anomalies
    df = leftjoin(df, baselines, on=:WellID)
    df = df[.!ismissing.(df.Baseline_m), :]
    df.Anomaly_m = compute_anomaly.(df.DepthToWater_m, df.Baseline_m)

    # Run QC
    df = run_bayesian_qc(df; cfg=cfg)

    # Filter to inliers for analysis
    df_clean = df[.!df.is_outlier, :]

    # Project coordinates
    xs, ys = project_to_albers(df_clean.Longitude, df_clean.Latitude; cfg=cfg)

    # Compute annual time series
    annual = compute_annual_timeseries(df_clean, xs, ys; cfg=cfg)

    # Get recent anomaly (IDW-weighted)
    recent_data = df_clean[df_clean.Year .>= cfg.recent_year, :]

    if nrow(recent_data) < 3
        # Fall back to last 5 years of data
        max_year = maximum(df_clean.Year)
        recent_data = df_clean[df_clean.Year .>= max_year - 5, :]
    end

    if nrow(recent_data) >= 3
        xs_recent, ys_recent = project_to_albers(recent_data.Longitude, recent_data.Latitude; cfg=cfg)
        anomaly_idw = idw_interpolate_single(xs_recent, ys_recent, recent_data.Anomaly_m;
                                              idp=cfg.idw_power, nmax=cfg.idw_neighbors)
        anomaly_simple = mean(recent_data.Anomaly_m)
        anomaly_sd = std(recent_data.Anomaly_m)
    else
        anomaly_idw = mean(df_clean.Anomaly_m)
        anomaly_simple = anomaly_idw
        anomaly_sd = std(df_clean.Anomaly_m)
    end

    # Estimate storage
    storage = estimate_storage(
        anomaly_idw, area_km2, sy_prior;
        aquifer_name=aquifer,
        anomaly_sd_m=anomaly_sd * 0.1,  # Reduce anomaly uncertainty
        n_samples=cfg.mc_samples
    )

    # Compute trend
    if nrow(annual) >= 3
        trend = compute_trend(annual.Year, annual.MeanAnomaly_m)
    else
        trend = (slope_decade=NaN, intercept=NaN, r2=NaN)
    end

    # Classify quality
    n_wells = length(unique(df_clean.WellID))
    n_years = length(unique(df_clean.Year))
    cv = abs(storage.storage_sd_km3 / storage.storage_mean_km3)
    quality = classify_quality(n_years, n_wells, trend.r2, cv)

    year_range = (minimum(df_clean.Year), maximum(df_clean.Year))

    return AquiferResult(
        aquifer,
        n_wells,
        n_years,
        year_range,
        trend.slope_decade,
        trend.r2,
        storage,
        quality,
        annual
    )
end

"""
    compute_annual_timeseries(df::DataFrame, xs::Vector, ys::Vector;
                              cfg::GWConfig=DEFAULT_CONFIG) -> DataFrame

Compute annual statistics with IDW weighting.
"""
function compute_annual_timeseries(
    df::DataFrame,
    xs::AbstractVector{<:Real},
    ys::AbstractVector{<:Real};
    cfg::GWConfig=DEFAULT_CONFIG
)
    years = sort(unique(df.Year))
    n_years = length(years)

    annual = DataFrame(
        Year = years,
        N_Wells = zeros(Int, n_years),
        MeanAnomaly_m = fill(NaN, n_years),
        MedianAnomaly_m = fill(NaN, n_years),
        IDWAnomaly_m = fill(NaN, n_years),
        SD_m = fill(NaN, n_years),
        Q05_m = fill(NaN, n_years),
        Q25_m = fill(NaN, n_years),
        Q75_m = fill(NaN, n_years),
        Q95_m = fill(NaN, n_years)
    )

    for (i, yr) in enumerate(years)
        idx = df.Year .== yr
        yr_data = df[idx, :]
        yr_x = xs[idx]
        yr_y = ys[idx]

        n = nrow(yr_data)
        annual.N_Wells[i] = n

        if n < 1
            continue
        end

        anom = yr_data.Anomaly_m

        annual.MeanAnomaly_m[i] = mean(anom)
        annual.MedianAnomaly_m[i] = median(anom)
        annual.SD_m[i] = n > 1 ? std(anom) : 0.0

        if n >= 5
            annual.Q05_m[i] = quantile(anom, 0.05)
            annual.Q25_m[i] = quantile(anom, 0.25)
            annual.Q75_m[i] = quantile(anom, 0.75)
            annual.Q95_m[i] = quantile(anom, 0.95)
        end

        # IDW mean
        if n >= 3
            annual.IDWAnomaly_m[i] = idw_interpolate_single(yr_x, yr_y, anom;
                                                            idp=cfg.idw_power,
                                                            nmax=min(cfg.idw_neighbors, n))
        else
            annual.IDWAnomaly_m[i] = mean(anom)
        end
    end

    return annual
end

"""
    run_storage_analysis(data::DataFrame; cfg::GWConfig=DEFAULT_CONFIG,
                         aquifer_areas::Dict{String, Float64}=Dict()) -> Vector{AquiferResult}

Run storage analysis for all aquifers in the data.

# Arguments
- `data`: Harmonized groundwater data
- `cfg`: Configuration parameters
- `aquifer_areas`: Dictionary of aquifer name => area (km²)

# Returns
Vector of AquiferResult for each aquifer
"""
function run_storage_analysis(
    data::DataFrame;
    cfg::GWConfig=DEFAULT_CONFIG,
    aquifer_areas::Dict{String, Float64}=Dict{String, Float64}()
)
    aquifers = get_aquifer_list(data)
    results = AquiferResult[]

    @info "Analyzing $(length(aquifers)) aquifers..."

    for aquifer in aquifers
        # Get area
        area = get(aquifer_areas, aquifer, NaN)
        if isnan(area)
            @warn "No area for $aquifer, skipping"
            continue
        end

        # Get aquifer data
        aq_data = get_wells_by_aquifer(data, aquifer)

        # Run analysis
        result = run_aquifer_analysis(aq_data, aquifer, area; cfg=cfg)

        if !isnothing(result)
            push!(results, result)
            s = result.storage
            @info "[$aquifer] $(result.n_years) yrs | ΔS=$(round(s.storage_mean_km3, digits=0))±$(round(s.storage_sd_km3, digits=0)) km³"
        end
    end

    return results
end

"""
    summarize_results(results::Vector{AquiferResult}) -> DataFrame

Create summary table from results.
"""
function summarize_results(results::Vector{AquiferResult})
    n = length(results)

    summary = DataFrame(
        Aquifer = [r.aquifer for r in results],
        Area_km2 = [r.storage.area_km2 for r in results],
        N_Wells = [r.n_wells for r in results],
        N_Years = [r.n_years for r in results],
        Year_Start = [r.year_range[1] for r in results],
        Year_End = [r.year_range[2] for r in results],
        Anomaly_m = [r.storage.anomaly_m for r in results],
        Sy_Mean = [r.storage.sy_prior.mean for r in results],
        Sy_SD = [r.storage.sy_prior.sd for r in results],
        Storage_km3 = [r.storage.storage_mean_km3 for r in results],
        Storage_SD_km3 = [r.storage.storage_sd_km3 for r in results],
        Storage_Q05_km3 = [r.storage.ci_95[1] for r in results],
        Storage_Q95_km3 = [r.storage.ci_95[2] for r in results],
        Trend_m_dec = [r.trend_m_decade for r in results],
        Trend_R2 = [r.trend_r2 for r in results],
        Quality = [String(r.quality_tier) for r in results],
        Confinement = [String(r.storage.sy_prior.confinement) for r in results]
    )

    # Sort by storage magnitude
    sort!(summary, :Storage_km3)

    return summary
end

"""
    export_results(results::Vector{AquiferResult}, output_dir::String)

Export results to CSV files.
"""
function export_results(results::Vector{AquiferResult}, output_dir::String)
    mkpath(output_dir)

    # Summary table
    summary = summarize_results(results)
    CSV.write(joinpath(output_dir, "storage_summary.csv"), summary)

    # Total CONUS
    total_mean = sum(summary.Storage_km3)
    total_sd = sqrt(sum(summary.Storage_SD_km3 .^ 2))
    ci_lo = total_mean - 1.96 * total_sd
    ci_hi = total_mean + 1.96 * total_sd

    totals = DataFrame(
        Metric = ["Total_km3", "SD_km3", "CI95_lo_km3", "CI95_hi_km3", "N_Aquifers"],
        Value = [total_mean, total_sd, ci_lo, ci_hi, length(results)]
    )
    CSV.write(joinpath(output_dir, "total_conus.csv"), totals)

    # Annual data per aquifer
    annual_dir = joinpath(output_dir, "annual")
    mkpath(annual_dir)

    for r in results
        CSV.write(joinpath(annual_dir, "$(r.aquifer)_annual.csv"), r.annual_data)
    end

    @info "Results exported to $output_dir"
    @info "  - storage_summary.csv"
    @info "  - total_conus.csv"
    @info "  - annual/<Aquifer>_annual.csv"

    return nothing
end

"""
    print_summary(results::Vector{AquiferResult})

Print summary to console.
"""
function print_summary(results::Vector{AquiferResult})
    summary = summarize_results(results)

    total_mean = sum(summary.Storage_km3)
    total_sd = sqrt(sum(summary.Storage_SD_km3 .^ 2))

    println("="^70)
    println("GROUNDWATER STORAGE ANALYSIS SUMMARY")
    println("="^70)
    println()
    println("Total aquifers: $(length(results))")
    println("Total CONUS storage change: $(round(total_mean, digits=0)) ± $(round(total_sd, digits=0)) km³")
    println("95% CI: [$(round(total_mean - 1.96*total_sd, digits=0)), $(round(total_mean + 1.96*total_sd, digits=0))] km³")
    println()
    println("Top 5 depleted aquifers:")
    for i in 1:min(5, nrow(summary))
        r = summary[i, :]
        println("  $(i). $(r.Aquifer): $(round(r.Storage_km3, digits=0)) ± $(round(r.Storage_SD_km3, digits=0)) km³")
    end
    println("="^70)
end

export print_summary

#!/usr/bin/env julia
# =============================================================================
# TREND BREAKPOINT ANALYSIS
# =============================================================================
#
# Detect acceleration/deceleration periods in groundwater depletion
# using piecewise linear regression.
#
# =============================================================================

using Pkg
Pkg.activate("/data2/GWProject")

using DataFrames, CSV, Statistics

const INPUT_DIR = "/data2/GWProject/outputs/BayesianSpatial"
const OUTPUT_DIR = "/data2/GWProject/outputs/TrendAnalysis"
mkpath(OUTPUT_DIR)

# =============================================================================
# Piecewise Linear Regression
# =============================================================================

function fit_piecewise_linear(years::Vector{Float64}, values::Vector{Float64}, breakpoint::Int)
    """Fit piecewise linear model with single breakpoint"""
    n = length(years)

    # Split data
    mask1 = years .<= breakpoint
    mask2 = years .> breakpoint

    n1, n2 = sum(mask1), sum(mask2)
    if n1 < 5 || n2 < 5
        return (slope1=NaN, slope2=NaN, r2=NaN, rss=Inf)
    end

    y1, v1 = years[mask1], values[mask1]
    y2, v2 = years[mask2], values[mask2]

    # Fit segment 1
    x1_mean, v1_mean = mean(y1), mean(v1)
    Sxx1 = sum((y1 .- x1_mean).^2)
    Sxy1 = sum((y1 .- x1_mean) .* (v1 .- v1_mean))
    slope1 = Sxy1 / Sxx1
    int1 = v1_mean - slope1 * x1_mean

    # Fit segment 2
    x2_mean, v2_mean = mean(y2), mean(v2)
    Sxx2 = sum((y2 .- x2_mean).^2)
    Sxy2 = sum((y2 .- x2_mean) .* (v2 .- v2_mean))
    slope2 = Sxy2 / Sxx2
    int2 = v2_mean - slope2 * x2_mean

    # Residual sum of squares
    pred1 = int1 .+ slope1 .* y1
    pred2 = int2 .+ slope2 .* y2
    rss = sum((v1 .- pred1).^2) + sum((v2 .- pred2).^2)

    # Total sum of squares
    v_mean = mean(values)
    tss = sum((values .- v_mean).^2)

    r2 = 1 - rss / tss

    return (slope1=slope1*10, slope2=slope2*10, r2=r2, rss=rss)
end

function find_optimal_breakpoint(years::Vector{Float64}, values::Vector{Float64})
    """Find breakpoint that minimizes RSS"""

    # Single linear fit first
    n = length(years)
    x_mean, y_mean = mean(years), mean(values)
    Sxx = sum((years .- x_mean).^2)
    Sxy = sum((years .- x_mean) .* (values .- y_mean))
    slope_single = Sxy / Sxx
    int_single = y_mean - slope_single * x_mean
    pred_single = int_single .+ slope_single .* years
    rss_single = sum((values .- pred_single).^2)
    tss = sum((values .- y_mean).^2)
    r2_single = 1 - rss_single / tss

    # Test breakpoints
    min_year = Int(minimum(years)) + 10
    max_year = Int(maximum(years)) - 10

    if max_year <= min_year
        return (breakpoint=0, slope1=slope_single*10, slope2=NaN,
                r2_piecewise=r2_single, r2_single=r2_single,
                improvement=0.0, acceleration=NaN)
    end

    best_bp = 0
    best_rss = Inf
    best_result = nothing

    for bp in min_year:5:max_year
        result = fit_piecewise_linear(years, values, bp)
        if result.rss < best_rss
            best_rss = result.rss
            best_bp = bp
            best_result = result
        end
    end

    if isnothing(best_result) || best_bp == 0
        return (breakpoint=0, slope1=slope_single*10, slope2=NaN,
                r2_piecewise=r2_single, r2_single=r2_single,
                improvement=0.0, acceleration=NaN)
    end

    # R² improvement
    improvement = best_result.r2 - r2_single

    # Acceleration (change in slope)
    acceleration = best_result.slope2 - best_result.slope1

    return (
        breakpoint = best_bp,
        slope1 = best_result.slope1,
        slope2 = best_result.slope2,
        r2_piecewise = best_result.r2,
        r2_single = r2_single,
        improvement = improvement,
        acceleration = acceleration
    )
end

# =============================================================================
# Main Analysis
# =============================================================================

println("=" ^ 70)
println("TREND BREAKPOINT ANALYSIS")
println("=" ^ 70)

# Load annual data for each aquifer
results = DataFrame()

aquifer_files = filter(f -> endswith(f, "_annual.csv"), readdir(INPUT_DIR))

for file in aquifer_files
    aquifer = replace(file, "_annual.csv" => "")

    df = CSV.read(joinpath(INPUT_DIR, file), DataFrame)

    # Get valid data
    valid = .!ismissing.(df.IDW_Mean) .& isfinite.(df.IDW_Mean)
    if sum(valid) < 20
        continue
    end

    years = Float64.(df.Year[valid])
    values = Float64.(df.IDW_Mean[valid])

    # Find optimal breakpoint
    bp_result = find_optimal_breakpoint(years, values)

    # Classify trend pattern
    if bp_result.improvement < 0.02
        pattern = "LINEAR"
    elseif bp_result.acceleration < -0.5
        pattern = "ACCELERATING_DEPLETION"
    elseif bp_result.acceleration > 0.5
        pattern = "DECELERATING_DEPLETION"
    elseif bp_result.slope1 < 0 && bp_result.slope2 > 0
        pattern = "RECOVERY"
    elseif bp_result.slope1 > 0 && bp_result.slope2 < 0
        pattern = "REVERSAL_TO_DEPLETION"
    else
        pattern = "STABLE"
    end

    push!(results, (
        Aquifer = aquifer,
        N_Years = sum(valid),
        Year_Start = Int(minimum(years)),
        Year_End = Int(maximum(years)),
        Breakpoint = bp_result.breakpoint,
        Slope_Before_m_dec = round(bp_result.slope1, digits=2),
        Slope_After_m_dec = round(bp_result.slope2, digits=2),
        Acceleration_m_dec2 = round(bp_result.acceleration, digits=2),
        R2_Single = round(bp_result.r2_single, digits=3),
        R2_Piecewise = round(bp_result.r2_piecewise, digits=3),
        R2_Improvement = round(bp_result.improvement, digits=3),
        Pattern = pattern
    ))

    if bp_result.improvement >= 0.02
        println("  $(rpad(aquifer, 35)) BP: $(bp_result.breakpoint) | $(round(bp_result.slope1, digits=1)) → $(round(bp_result.slope2, digits=1)) m/dec | $pattern")
    end
end

# =============================================================================
# Summary
# =============================================================================

println("\n" * "=" ^ 70)
println("BREAKPOINT SUMMARY")
println("=" ^ 70)

# Count patterns
pattern_counts = combine(groupby(results, :Pattern), nrow => :Count)
sort!(pattern_counts, :Count, rev=true)

println("\nTrend Patterns:")
for row in eachrow(pattern_counts)
    println("  $(rpad(row.Pattern, 25)) $(row.Count) aquifers")
end

# Significant breakpoints (R² improvement > 5%)
sig_bp = results[results.R2_Improvement .>= 0.05, :]
sort!(sig_bp, :R2_Improvement, rev=true)

println("\n" * "-" ^ 70)
println("SIGNIFICANT BREAKPOINTS (R² improvement ≥ 5%)")
println("-" ^ 70)

for row in eachrow(sig_bp)
    println("  $(rpad(row.Aquifer, 30)) $(row.Breakpoint): $(row.Slope_Before_m_dec) → $(row.Slope_After_m_dec) m/dec (ΔR²=$(row.R2_Improvement))")
end

# Accelerating depletion
accel = results[results.Pattern .== "ACCELERATING_DEPLETION", :]
if nrow(accel) > 0
    println("\n" * "-" ^ 70)
    println("ACCELERATING DEPLETION (rates increasing)")
    println("-" ^ 70)
    for row in eachrow(accel)
        println("  $(row.Aquifer): $(row.Slope_Before_m_dec) → $(row.Slope_After_m_dec) m/dec after $(row.Breakpoint)")
    end
end

# Recovery
recovery = results[results.Pattern .== "RECOVERY", :]
if nrow(recovery) > 0
    println("\n" * "-" ^ 70)
    println("RECOVERY (depletion → recovery)")
    println("-" ^ 70)
    for row in eachrow(recovery)
        println("  $(row.Aquifer): $(row.Slope_Before_m_dec) → $(row.Slope_After_m_dec) m/dec after $(row.Breakpoint)")
    end
end

# =============================================================================
# Save Results
# =============================================================================

sort!(results, :Acceleration_m_dec2)
CSV.write(joinpath(OUTPUT_DIR, "trend_breakpoints.csv"), results)

println("\n" * "=" ^ 70)
println("OUTPUT: $(joinpath(OUTPUT_DIR, "trend_breakpoints.csv"))")
println("=" ^ 70)

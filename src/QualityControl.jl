# =============================================================================
# Quality Control Methods
# =============================================================================

"""
    WellClass

Classification of well behavior based on residual analysis.
"""
@enum WellClass begin
    INSUFFICIENT_DATA
    LOW_QUALITY
    LIKELY_CONFINED
    LIKELY_UNCONFINED
    INDETERMINATE
end

"""
    detect_outliers(values::AbstractVector; z_threshold::Real=3.5) -> BitVector

Detect outliers using MAD-based Z-score.

# Arguments
- `values`: Vector of values to check
- `z_threshold`: Z-score threshold (default: 3.5)

# Returns
BitVector where true = outlier
"""
function detect_outliers(values::AbstractVector{<:Real}; z_threshold::Real=3.5)
    valid = isfinite.(values)
    v = values[valid]

    if length(v) < 3
        return falses(length(values))
    end

    # MAD-based robust Z-score
    med = median(v)
    mad = median(abs.(v .- med)) * 1.4826  # Consistency constant for normal

    if mad < 1e-10
        return falses(length(values))
    end

    z_scores = abs.(values .- med) ./ mad
    z_scores[.!valid] .= 0

    return z_scores .> z_threshold
end

"""
    robust_mad(x::AbstractVector; constant::Real=1.4826) -> Float64

Compute Median Absolute Deviation with consistency constant.
"""
function robust_mad(x::AbstractVector{<:Real}; constant::Real=1.4826)
    v = filter(isfinite, x)
    if isempty(v)
        return NaN
    end
    med = median(v)
    return median(abs.(v .- med)) * constant
end

"""
    t_mixture_posterior(residual::Real, sigma::Real, pi0::Real;
                        df_in::Real=7.0, df_out::Real=1.0, scale_mult::Real=3.0) -> Float64

Compute posterior probability of being an outlier using t-mixture model.

# Model
- Inlier: t(df_in) with scale sigma
- Outlier: t(df_out) with scale sigma × scale_mult

# Returns
P(outlier | residual)
"""
function t_mixture_posterior(
    residual::Real,
    sigma::Real,
    pi0::Real;
    df_in::Real=7.0,
    df_out::Real=1.0,
    scale_mult::Real=3.0
)
    if !isfinite(residual) || !isfinite(sigma) || sigma <= 0
        return 0.0
    end

    z_in = residual / sigma
    z_out = z_in / scale_mult

    # t-distribution densities
    dens_in = pdf(TDist(df_in), z_in)
    dens_out = pdf(TDist(df_out), z_out) / scale_mult

    # Posterior
    numer = pi0 * dens_out
    denom = numer + (1.0 - pi0) * dens_in

    return denom > 0 ? numer / denom : 0.0
end

"""
    run_bayesian_qc(df::DataFrame; cfg::GWConfig=DEFAULT_CONFIG) -> DataFrame

Run Bayesian quality control on groundwater data.

# Process
1. Fit robust trend per well
2. Compute residuals
3. Calculate outlier posterior probabilities
4. Flag outliers and compute QC weights

# Returns
DataFrame with added columns:
- fitted, residual, p_outlier, is_outlier, w_qc
"""
function run_bayesian_qc(df::DataFrame; cfg::GWConfig=DEFAULT_CONFIG)
    df = copy(df)

    # Initialize columns
    df.fitted = fill(NaN, nrow(df))
    df.residual = fill(NaN, nrow(df))
    df.p_outlier = fill(0.0, nrow(df))
    df.is_outlier = falses(nrow(df))
    df.w_qc = fill(1.0, nrow(df))

    # Get source-specific prior outlier probabilities
    source_pi0 = Dict(
        "NWIS" => 0.005,
        "CA" => 0.012,
        "TX" => 0.015
    )
    default_pi0 = cfg.outlier_prior

    # Process each well
    wells = unique(df.WellID)

    for well in wells
        idx = findall(df.WellID .== well)
        well_data = df[idx, :]

        if nrow(well_data) < 3
            continue
        end

        # Get source prior
        src = first(well_data.Source)
        pi0 = get(source_pi0, src, default_pi0)

        # Fit OLS trend
        years = Float64.(well_data.Year)
        depths = well_data.DepthToWater_m

        year_mean = mean(years)
        depth_mean = mean(depths)
        Stt = sum((years .- year_mean).^2)
        Std = sum((years .- year_mean) .* (depths .- depth_mean))

        slope = Stt > 0 ? Std / Stt : 0.0
        intercept = depth_mean - slope * year_mean

        fitted = intercept .+ slope .* years
        resid = depths .- fitted

        # Robust scale
        sigma = robust_mad(resid)
        if !isfinite(sigma) || sigma < 0.05
            sigma = 0.05  # Floor
        end

        # Compute outlier posteriors
        p_out = [t_mixture_posterior(r, sigma, pi0) for r in resid]

        # Store results
        df.fitted[idx] = fitted
        df.residual[idx] = resid
        df.p_outlier[idx] = p_out
        df.is_outlier[idx] = p_out .> 0.5
        df.w_qc[idx] = 1.0 .- p_out
    end

    n_outliers = sum(df.is_outlier)
    pct = round(100 * n_outliers / nrow(df), digits=2)
    @info "Detected $n_outliers outliers ($pct%)"

    return df
end

"""
    classify_well(residuals::AbstractVector, years::AbstractVector;
                  amp_small::Real=1.0, amp_large::Real=3.0,
                  r2_min::Real=0.70, n_min::Int=3) -> WellClass

Classify well behavior based on residual patterns.

# Criteria
- LIKELY_CONFINED: Low amplitude (P90-P10 < amp_small), smooth, high R²
- LIKELY_UNCONFINED: High amplitude (P90-P10 > amp_large)
- INDETERMINATE: Everything else
"""
function classify_well(
    residuals::AbstractVector{<:Real},
    years::AbstractVector{<:Real};
    amp_small::Real=1.0,
    amp_large::Real=3.0,
    r2_min::Real=0.70,
    n_min::Int=3
)
    # Filter valid
    valid = isfinite.(residuals) .& isfinite.(years)
    r = residuals[valid]
    t = years[valid]

    n_years = length(unique(t))

    if n_years < n_min
        return INSUFFICIENT_DATA
    end

    # Amplitude: P90 - P10
    amp = quantile(r, 0.90) - quantile(r, 0.10)

    if !isfinite(amp)
        return INDETERMINATE
    end

    if amp >= amp_large
        return LIKELY_UNCONFINED
    end

    # For confined classification, need smooth residuals and good trend fit
    if length(r) >= 2
        ord = sortperm(t)
        dr = diff(r[ord])
        dt = diff(t[ord])
        dt[dt .== 0] .= 1
        med_dy = median(abs.(dr ./ dt))
    else
        med_dy = NaN
    end

    # R² would need original values - simplified here
    if amp <= amp_small && isfinite(med_dy) && med_dy <= 0.25
        return LIKELY_CONFINED
    end

    return INDETERMINATE
end

"""
    winsorize(x::AbstractVector, prob::Real=0.05) -> Vector{Float64}

Winsorize extreme values.

# Arguments
- `x`: Input vector
- `prob`: Probability for each tail (default: 0.05 = 5% each tail)
"""
function winsorize(x::AbstractVector{<:Real}, prob::Real=0.05)
    lo = quantile(filter(isfinite, x), prob)
    hi = quantile(filter(isfinite, x), 1 - prob)
    return clamp.(x, lo, hi)
end

export robust_mad, t_mixture_posterior, winsorize

# =============================================================================
# Spatial Methods: Projection, IDW, Grid Operations
# =============================================================================

"""
    project_to_albers(lon::Real, lat::Real; cfg::GWConfig=DEFAULT_CONFIG) -> Tuple{Float64, Float64}

Project WGS84 coordinates to NAD83/Conus Albers Equal-Area (EPSG:5070).

# Arguments
- `lon`: Longitude in degrees (WGS84)
- `lat`: Latitude in degrees (WGS84)
- `cfg`: Configuration with projection parameters

# Returns
Tuple (x, y) in meters (Albers coordinates)

# Example
```julia
x, y = project_to_albers(-96.0, 35.0)
```
"""
function project_to_albers(lon::Real, lat::Real; cfg::GWConfig=DEFAULT_CONFIG)
    # NAD83/Conus Albers parameters (EPSG:5070)
    lat0 = cfg.proj_lat0   # 23.0
    lon0 = cfg.proj_lon0   # -96.0
    lat1 = cfg.proj_lat1   # 29.5
    lat2 = cfg.proj_lat2   # 45.5

    # Convert to radians
    φ = deg2rad(lat)
    λ = deg2rad(lon)
    φ0 = deg2rad(lat0)
    λ0 = deg2rad(lon0)
    φ1 = deg2rad(lat1)
    φ2 = deg2rad(lat2)

    # GRS80 ellipsoid
    a = 6378137.0
    e2 = 0.00669438
    e = sqrt(e2)

    # Helper functions
    q(φ) = (1 - e2) * (sin(φ) / (1 - e2 * sin(φ)^2) -
           (1 / (2e)) * log((1 - e * sin(φ)) / (1 + e * sin(φ))))
    m(φ) = cos(φ) / sqrt(1 - e2 * sin(φ)^2)

    q0 = q(φ0)
    q1 = q(φ1)
    q2 = q(φ2)
    m1 = m(φ1)
    m2 = m(φ2)

    n = (m1^2 - m2^2) / (q2 - q1)
    C = m1^2 + n * q1
    ρ0 = a * sqrt(C - n * q0) / n

    qφ = q(φ)
    ρ = a * sqrt(C - n * qφ) / n
    θ = n * (λ - λ0)

    x = ρ * sin(θ)
    y = ρ0 - ρ * cos(θ)

    return (x, y)
end

"""
    project_to_albers(lons::AbstractVector, lats::AbstractVector; cfg::GWConfig=DEFAULT_CONFIG) -> Tuple{Vector{Float64}, Vector{Float64}}

Vectorized projection for arrays of coordinates.
"""
function project_to_albers(lons::AbstractVector, lats::AbstractVector; cfg::GWConfig=DEFAULT_CONFIG)
    n = length(lons)
    @assert length(lats) == n "lon and lat must have same length"

    xs = Vector{Float64}(undef, n)
    ys = Vector{Float64}(undef, n)

    for i in 1:n
        xs[i], ys[i] = project_to_albers(lons[i], lats[i]; cfg=cfg)
    end

    return (xs, ys)
end

"""
    inverse_project(x::Real, y::Real; cfg::GWConfig=DEFAULT_CONFIG) -> Tuple{Float64, Float64}

Inverse projection from Albers to WGS84 (approximate).

# Returns
Tuple (lon, lat) in degrees
"""
function inverse_project(x::Real, y::Real; cfg::GWConfig=DEFAULT_CONFIG)
    # Iterative inverse projection (simplified)
    lat0 = cfg.proj_lat0
    lon0 = cfg.proj_lon0
    lat1 = cfg.proj_lat1
    lat2 = cfg.proj_lat2

    φ0 = deg2rad(lat0)
    λ0 = deg2rad(lon0)
    φ1 = deg2rad(lat1)
    φ2 = deg2rad(lat2)

    a = 6378137.0
    e2 = 0.00669438
    e = sqrt(e2)

    q(φ) = (1 - e2) * (sin(φ) / (1 - e2 * sin(φ)^2) -
           (1 / (2e)) * log((1 - e * sin(φ)) / (1 + e * sin(φ))))
    m(φ) = cos(φ) / sqrt(1 - e2 * sin(φ)^2)

    q0 = q(φ0)
    q1 = q(φ1)
    q2 = q(φ2)
    m1 = m(φ1)
    m2 = m(φ2)

    n = (m1^2 - m2^2) / (q2 - q1)
    C = m1^2 + n * q1
    ρ0 = a * sqrt(C - n * q0) / n

    ρ = sqrt(x^2 + (ρ0 - y)^2)
    θ = atan(x, ρ0 - y)

    λ = λ0 + θ / n
    qφ = (C - (ρ * n / a)^2) / n

    # Iterative solve for φ
    φ = asin(qφ / 2)
    for _ in 1:10
        sinφ = sin(φ)
        φ_new = φ + (1 - e2 * sinφ^2)^2 / (2 * cos(φ)) *
                (qφ / (1 - e2) - sinφ / (1 - e2 * sinφ^2) +
                 (1 / (2e)) * log((1 - e * sinφ) / (1 + e * sinφ)))
        if abs(φ_new - φ) < 1e-12
            break
        end
        φ = φ_new
    end

    return (rad2deg(λ), rad2deg(φ))
end

"""
    idw_interpolate(x_obs, y_obs, z_obs, x_pred, y_pred;
                    idp=2, nmax=12, min_dist=1.0) -> Vector{Float64}

Inverse Distance Weighting interpolation.

# Arguments
- `x_obs, y_obs`: Observation coordinates (projected, meters)
- `z_obs`: Observation values
- `x_pred, y_pred`: Prediction coordinates
- `idp`: Power parameter (default: 2)
- `nmax`: Maximum number of neighbors (default: 12)
- `min_dist`: Minimum distance to avoid singularities (default: 1.0 m)

# Returns
Vector of interpolated values at prediction points

# Example
```julia
z_pred = idw_interpolate(x_wells, y_wells, anomaly, x_grid, y_grid; idp=2, nmax=12)
```
"""
function idw_interpolate(
    x_obs::AbstractVector{<:Real},
    y_obs::AbstractVector{<:Real},
    z_obs::AbstractVector{<:Real},
    x_pred::AbstractVector{<:Real},
    y_pred::AbstractVector{<:Real};
    idp::Real=2,
    nmax::Int=12,
    min_dist::Real=1.0
)
    n_obs = length(x_obs)
    n_pred = length(x_pred)

    @assert length(y_obs) == n_obs "x_obs and y_obs must have same length"
    @assert length(z_obs) == n_obs "z_obs must have same length as coordinates"
    @assert length(y_pred) == n_pred "x_pred and y_pred must have same length"

    # Build KD-tree for fast neighbor search
    coords = hcat(x_obs, y_obs)'
    tree = KDTree(coords)

    # Number of neighbors to use
    k = min(nmax, n_obs)

    z_pred = Vector{Float64}(undef, n_pred)

    for i in 1:n_pred
        query = [x_pred[i], y_pred[i]]
        idxs, dists = knn(tree, query, k, true)  # sortres=true

        # Apply minimum distance
        dists = max.(dists, min_dist)

        # Compute weights
        weights = 1.0 ./ (dists .^ idp)

        # Weighted average
        z_vals = z_obs[idxs]
        z_pred[i] = sum(weights .* z_vals) / sum(weights)
    end

    return z_pred
end

"""
    idw_interpolate_single(x_obs, y_obs, z_obs; idp=2, nmax=12) -> Float64

Compute IDW mean of observations (interpolate at centroid).
"""
function idw_interpolate_single(
    x_obs::AbstractVector{<:Real},
    y_obs::AbstractVector{<:Real},
    z_obs::AbstractVector{<:Real};
    idp::Real=2,
    nmax::Int=12
)
    # Interpolate at centroid
    x_cent = mean(x_obs)
    y_cent = mean(y_obs)

    result = idw_interpolate(x_obs, y_obs, z_obs, [x_cent], [y_cent]; idp=idp, nmax=nmax)
    return result[1]
end

"""
    compute_nyquist_spacing(x::AbstractVector, y::AbstractVector) -> Float64

Compute Nyquist-optimal grid spacing based on nearest neighbor distances.

# Returns
Recommended grid spacing = median(NN_distance) / 2
"""
function compute_nyquist_spacing(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(x)
    if n < 2
        return NaN
    end

    coords = hcat(x, y)'
    tree = KDTree(coords)

    # Find nearest neighbor for each point
    nn_dists = Float64[]
    for i in 1:n
        query = [x[i], y[i]]
        idxs, dists = knn(tree, query, 2, true)
        # dists[1] is 0 (self), dists[2] is nearest neighbor
        if length(dists) >= 2
            push!(nn_dists, dists[2])
        end
    end

    # Nyquist spacing
    return median(nn_dists) / 2
end

"""
    assign_wells_to_grid(lons::AbstractVector, lats::AbstractVector,
                         resolution::Real) -> Tuple{Vector{Int}, Vector{Int}}

Assign wells to grid cells.

# Arguments
- `lons, lats`: Well coordinates in degrees
- `resolution`: Grid resolution in degrees

# Returns
Tuple of (row_indices, col_indices)
"""
function assign_wells_to_grid(
    lons::AbstractVector{<:Real},
    lats::AbstractVector{<:Real},
    resolution::Real
)
    # Grid indices (0-based then convert to 1-based)
    cols = floor.(Int, (lons .+ 180) ./ resolution) .+ 1
    rows = floor.(Int, (lats .+ 90) ./ resolution) .+ 1

    return (rows, cols)
end

"""
    compute_pairwise_distances(x::AbstractVector, y::AbstractVector) -> Matrix{Float64}

Compute pairwise Euclidean distances.
"""
function compute_pairwise_distances(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(x)
    D = Matrix{Float64}(undef, n, n)

    for i in 1:n
        for j in 1:n
            D[i, j] = sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)
        end
    end

    return D
end

export idw_interpolate_single, compute_pairwise_distances

# =============================================================================
# Data Input/Output Functions
# =============================================================================

"""
    load_harmonized_data(filepath::String; cfg::GWConfig=DEFAULT_CONFIG) -> DataFrame

Load harmonized groundwater data from CSV file.

# Arguments
- `filepath`: Path to harmonized CSV file
- `cfg`: Configuration parameters

# Returns
DataFrame with standardized columns:
- WellID, Latitude, Longitude, Year, Season
- DepthToWater_m, Source, Aquifer

# Example
```julia
df = load_harmonized_data("/data/GW_CONUS_Harmonized.csv")
```
"""
function load_harmonized_data(filepath::String; cfg::GWConfig=DEFAULT_CONFIG)
    @info "Loading harmonized data from $filepath"

    df = CSV.read(filepath, DataFrame; missingstring=["NA", "NaN", ""])

    # Standardize column names
    if :DepthToWater_ft in propertynames(df) && !(:DepthToWater_m in propertynames(df))
        df.DepthToWater_m = df.DepthToWater_ft .* 0.3048
    end

    # Extract year if needed
    if !(:Year in propertynames(df)) && :MeasurementDate in propertynames(df)
        df.Year = year.(df.MeasurementDate)
    end

    # Filter by bounds (CONUS)
    valid_mask = (df.Latitude .>= 24.0) .& (df.Latitude .<= 50.0) .&
                 (df.Longitude .>= -125.0) .& (df.Longitude .<= -66.0)
    df = df[valid_mask, :]

    # Filter by year range
    df = df[(df.Year .>= cfg.year_min) .& (df.Year .<= cfg.year_max), :]

    # Remove invalid depths
    df = df[df.DepthToWater_m .> 0, :]
    df = df[df.DepthToWater_m .< 1000, :]  # Max reasonable depth

    @info "Loaded $(nrow(df)) records from $(length(unique(df.WellID))) wells"

    return df
end

"""
    filter_by_season(df::DataFrame, seasons::Vector{String}) -> DataFrame

Filter data to specified seasons.

# Arguments
- `df`: Input DataFrame with Season column
- `seasons`: Vector of season names (e.g., ["Winter", "Spring"])
"""
function filter_by_season(df::DataFrame, seasons::Vector{String})
    if :Season in propertynames(df)
        return df[df.Season .∈ Ref(seasons), :]
    else
        @warn "No Season column found, returning all data"
        return df
    end
end

"""
    filter_by_years(df::DataFrame, year_start::Int, year_end::Int) -> DataFrame

Filter data to specified year range.
"""
function filter_by_years(df::DataFrame, year_start::Int, year_end::Int)
    return df[(df.Year .>= year_start) .& (df.Year .<= year_end), :]
end

"""
    compute_well_statistics(df::DataFrame) -> DataFrame

Compute per-well statistics.

# Returns
DataFrame with one row per well containing:
- WellID, Latitude, Longitude, Aquifer, Source
- n_obs, n_years, year_min, year_max
- depth_mean, depth_median, depth_sd, depth_iqr
"""
function compute_well_statistics(df::DataFrame)
    gdf = groupby(df, :WellID)

    stats = combine(gdf,
        :Latitude => median => :Latitude,
        :Longitude => median => :Longitude,
        :Aquifer => first => :Aquifer,
        :Source => first => :Source,
        nrow => :n_obs,
        :Year => (y -> length(unique(y))) => :n_years,
        :Year => minimum => :year_min,
        :Year => maximum => :year_max,
        :DepthToWater_m => mean => :depth_mean,
        :DepthToWater_m => median => :depth_median,
        :DepthToWater_m => std => :depth_sd,
        :DepthToWater_m => iqr => :depth_iqr
    )

    return stats
end

"""
    load_aquifer_boundaries(gpkg_path::String) -> DataFrame

Load EPA Principal Aquifer boundaries from GeoPackage.

Note: Requires GeoDataFrames.jl to be installed separately.
For most use cases, use pre-computed aquifer areas instead.

# Returns
DataFrame with aquifer names and areas.
"""
function load_aquifer_boundaries(gpkg_path::String)
    @warn "GeoDataFrames not loaded. Use pre-computed aquifer areas dictionary instead."
    @warn "See examples/run_analysis.jl for aquifer_areas dictionary."
    return DataFrame(Aquifer=String[], Area_km2=Float64[])
end

"""
    get_wells_by_aquifer(df::DataFrame, aquifer::String) -> DataFrame

Extract data for a specific aquifer.
"""
function get_wells_by_aquifer(df::DataFrame, aquifer::String)
    return df[df.Aquifer .== aquifer, :]
end

"""
    get_aquifer_list(df::DataFrame) -> Vector{String}

Get list of unique aquifers in the data.
"""
function get_aquifer_list(df::DataFrame)
    return sort(unique(df.Aquifer))
end

"""
    export_to_csv(df::DataFrame, filepath::String)

Export DataFrame to CSV with standard formatting.
"""
function export_to_csv(df::DataFrame, filepath::String)
    CSV.write(filepath, df)
    @info "Exported to $filepath"
end

export get_wells_by_aquifer, get_aquifer_list, export_to_csv

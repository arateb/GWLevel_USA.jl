# =============================================================================
# Configuration Types and Defaults
# =============================================================================

"""
    GWConfig

Configuration parameters for groundwater analysis.

# Fields
- `baseline_start::Int`: Start year for baseline period (default: 1940)
- `baseline_end::Int`: End year for baseline period (default: 1955)
- `valid_seasons::Vector{String}`: Seasons to include (default: ["Winter", "Spring"])
- `year_min::Int`: Minimum year to include (default: 1940)
- `year_max::Int`: Maximum year to include (default: 2025)
- `min_obs_per_well::Int`: Minimum observations per well (default: 5)
- `min_years_baseline::Int`: Minimum years in baseline period (default: 2)
- `min_wells_per_aquifer::Int`: Minimum wells per aquifer (default: 5)
- `outlier_zscore::Float64`: Z-score threshold for outliers (default: 3.5)
- `idw_power::Float64`: IDW interpolation power (default: 2.0)
- `idw_neighbors::Int`: Maximum neighbors for IDW (default: 12)
- `mc_samples::Int`: Monte Carlo samples (default: 10000)
- `grid_resolution::Float64`: Grid resolution in degrees (default: 0.25)
"""
Base.@kwdef struct GWConfig
    # Temporal parameters
    baseline_start::Int = 1940
    baseline_end::Int = 1955
    valid_seasons::Vector{String} = ["Winter", "Spring"]
    year_min::Int = 1940
    year_max::Int = 2025
    recent_year::Int = 2020

    # Data quality thresholds
    min_obs_per_well::Int = 5
    min_years_baseline::Int = 2
    min_wells_per_aquifer::Int = 5
    min_records_per_aquifer::Int = 100

    # Outlier detection
    outlier_zscore::Float64 = 3.5
    outlier_prior::Float64 = 0.01

    # Spatial parameters
    idw_power::Float64 = 2.0
    idw_neighbors::Int = 12
    grid_resolution::Float64 = 0.25

    # Monte Carlo
    mc_samples::Int = 10000

    # Projection (EPSG:5070 NAD83/Conus Albers)
    proj_lat0::Float64 = 23.0
    proj_lon0::Float64 = -96.0
    proj_lat1::Float64 = 29.5
    proj_lat2::Float64 = 45.5
end

"""
    SyPrior

Specific yield prior distribution for an aquifer.

# Fields
- `mean::Float64`: Prior mean
- `sd::Float64`: Prior standard deviation
- `reference::String`: Literature reference
- `confinement::Symbol`: :unconfined, :confined, or :mixed
"""
struct SyPrior
    mean::Float64
    sd::Float64
    reference::String
    confinement::Symbol

    function SyPrior(mean, sd, reference, confinement)
        @assert mean > 0 "Sy mean must be positive"
        @assert sd > 0 "Sy SD must be positive"
        @assert confinement in [:unconfined, :confined, :mixed] "Invalid confinement type"
        new(mean, sd, reference, confinement)
    end
end

"""
    AquiferInfo

Information about an aquifer including area and Sy prior.

# Fields
- `name::String`: Aquifer name (ID)
- `full_name::String`: Full aquifer name
- `area_km2::Float64`: Area in km²
- `sy_prior::SyPrior`: Specific yield prior
"""
struct AquiferInfo
    name::String
    full_name::String
    area_km2::Float64
    sy_prior::SyPrior
end

# Default configuration
const DEFAULT_CONFIG = GWConfig()

# Literature-based Sy priors
const DEFAULT_SY_PRIORS = Dict{String, SyPrior}(
    # High Plains
    "N_High_Plains" => SyPrior(0.15, 0.04, "Gutentag et al. (1984)", :unconfined),
    "C_S_High_Plains" => SyPrior(0.15, 0.04, "Gutentag et al. (1984)", :unconfined),
    "High_Plains" => SyPrior(0.15, 0.04, "McGuire (2017)", :unconfined),

    # Central Valley
    "Central_Valley" => SyPrior(0.12, 0.04, "Faunt (2009)", :unconfined),
    "Sacramento" => SyPrior(0.12, 0.04, "Faunt (2009)", :unconfined),
    "San_Joaquin" => SyPrior(0.12, 0.04, "Faunt (2009)", :unconfined),

    # Basin and Range
    "Basin_and_Range_Basin_Fill" => SyPrior(0.18, 0.04, "Pool & Coes (1999)", :unconfined),
    "Arizona_Alluvials" => SyPrior(0.18, 0.05, "Pool & Coes (1999)", :unconfined),

    # Mississippi System
    "Mississippi_River_Valley_Alluvial" => SyPrior(0.22, 0.05, "Clark & Hart (2009)", :unconfined),
    "Mississippi_Embayment" => SyPrior(0.12, 0.04, "Clark & Hart (2009)", :mixed),

    # Pacific Northwest
    "Snake" => SyPrior(0.08, 0.04, "Whitehead (1992)", :mixed),
    "Columbia_Plateau_Basin_fill" => SyPrior(0.12, 0.05, "Kahle et al. (2011)", :unconfined),
    "Columbia_Plateau_Basaltic_rock" => SyPrior(0.05, 0.03, "Kahle et al. (2011)", :mixed),

    # Gulf Coast
    "Coastal_Lowland" => SyPrior(0.08, 0.04, "Kasmarek et al. (2016)", :mixed),
    "Texas_Gulf_Coast" => SyPrior(0.06, 0.03, "Kasmarek et al. (2016)", :confined),

    # Atlantic Coastal Plain
    "N_Atlantic_Coastal_Plain" => SyPrior(0.12, 0.04, "Masterson et al. (2016)", :mixed),
    "SE_Coastal_Plain" => SyPrior(0.10, 0.04, "Campbell & Coes (2010)", :mixed),

    # Floridan (confined karst)
    "Floridan" => SyPrior(0.001, 0.0005, "Miller (1986)", :confined),

    # Glacial
    "NE_Glaciers" => SyPrior(0.18, 0.06, "Starn & Brown (2007)", :unconfined),
    "Great_lakes" => SyPrior(0.15, 0.05, "Grannemann et al. (2000)", :unconfined),

    # Bedrock (storativity for confined)
    "Cambrian_Ordovician" => SyPrior(0.0005, 0.0003, "Young (1992)", :confined),
    "Piedmont_and_Blue_Ridge" => SyPrior(0.01, 0.005, "Daniel (1989)", :mixed),
    "Ozark_Plateau" => SyPrior(0.02, 0.01, "Imes & Emmett (1994)", :mixed),

    # Texas
    "Carrizo_Wilcox" => SyPrior(0.0003, 0.0002, "Mace et al. (2000)", :confined),
    "Trinity" => SyPrior(0.0005, 0.0003, "Kelley et al. (2004)", :confined),
    "Edwards_Plateau" => SyPrior(0.035, 0.015, "Maclay & Small (1986)", :mixed),
    "Rio_Grande" => SyPrior(0.15, 0.05, "Heywood & Yager (2003)", :unconfined),

    # Other
    "U_Colorado" => SyPrior(0.08, 0.04, "Robson & Banta (1995)", :mixed),
    "N_Great_Plains" => SyPrior(0.12, 0.05, "Downey (1986)", :mixed),
    "Pennsylvania" => SyPrior(0.05, 0.03, "Default - fractured rock", :mixed),
)

"""
    get_sy_prior(aquifer_name::String) -> SyPrior

Get specific yield prior for an aquifer. Returns default if not found.
"""
function get_sy_prior(aquifer_name::String)
    if haskey(DEFAULT_SY_PRIORS, aquifer_name)
        return DEFAULT_SY_PRIORS[aquifer_name]
    else
        # Default for unknown aquifers
        return SyPrior(0.10, 0.05, "Default estimate", :mixed)
    end
end

export get_sy_prior

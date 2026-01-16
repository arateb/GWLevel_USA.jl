"""
    GWLevel_USA

Bayesian framework for estimating groundwater storage changes across the
conterminous United States using in-situ well observations.

# Main Features
- Load and preprocess harmonized groundwater data
- Assign wells to EPA Principal Aquifers
- Compute baseline depths and anomalies
- IDW spatial interpolation to correct monitoring bias
- Bayesian storage estimation with literature-based Sy priors
- Quality control with t-mixture outlier detection
- Monte Carlo uncertainty quantification

# Quick Start
```julia
using GWLevel_USA

# Load data
data = load_harmonized_data("/path/to/GW_CONUS_Harmonized.csv")

# Run full analysis
results = run_storage_analysis(data; baseline_period=1940:1955)

# Get summary
summary(results)
```

# Author
Dr. Ashraf Rateb, Bureau of Economic Geology, University of Texas at Austin
"""
module GWLevel_USA

using CSV
using DataFrames
using Dates
using Distributions
using LinearAlgebra
using NearestNeighbors
using Printf
using Random
using Statistics
using StatsBase

# Note: GeoDataFrames is optional for loading .gpkg files
# If needed, add to Project.toml: GeoDataFrames = "62cb38b5-d8d2-4862-a48e-6a340996859f"

# Include submodules
include("Config.jl")
include("DataIO.jl")
include("Spatial.jl")
include("Bayesian.jl")
include("QualityControl.jl")
include("Analysis.jl")

# Re-export from Config
export GWConfig, SyPrior, AquiferInfo
export DEFAULT_CONFIG, DEFAULT_SY_PRIORS

# Re-export from DataIO
export load_harmonized_data, load_aquifer_boundaries
export filter_by_season, filter_by_years
export compute_well_statistics

# Re-export from Spatial
export project_to_albers, inverse_project
export idw_interpolate, compute_nyquist_spacing
export assign_wells_to_grid

# Re-export from Bayesian
export compute_baseline, compute_anomaly
export estimate_storage, run_monte_carlo
export StorageResult, AquiferResult

# Re-export from QualityControl
export detect_outliers, run_bayesian_qc
export classify_well, WellClass

# Re-export from Analysis
export run_storage_analysis, run_aquifer_analysis
export summarize_results, export_results

end # module

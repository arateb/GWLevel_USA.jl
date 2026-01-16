#!/usr/bin/env julia
# =============================================================================
# Example: Run Full CONUS Groundwater Storage Analysis
# =============================================================================
#
# This script demonstrates how to use the GWLevel_USA package to analyze
# groundwater storage changes across the conterminous United States.
#
# Author: Dr. Ashraf Rateb, Bureau of Economic Geology, UT Austin
# =============================================================================

using Pkg
Pkg.activate(dirname(@__DIR__))

using GWLevel_USA
using CSV
using DataFrames

# =============================================================================
# Configuration
# =============================================================================

# Data paths - adjust to your system
DATA_DIR = "/data2/GWProject/data"
OUTPUT_DIR = "/data2/GWProject/outputs/GWLevel_USA"

HARMONIZED_CSV = joinpath(DATA_DIR, "CONUS/GW_CONUS_Harmonized_20251208.csv")
AQUIFER_GPKG = joinpath(DATA_DIR, "EPA_Aquifers_Oct2025.gpkg")

# Analysis configuration
cfg = GWConfig(
    baseline_start = 1940,
    baseline_end = 1955,
    valid_seasons = ["Winter", "Spring"],
    min_obs_per_well = 5,
    min_wells_per_aquifer = 5,
    idw_power = 2.0,
    idw_neighbors = 12,
    mc_samples = 10000
)

# =============================================================================
# Load Data
# =============================================================================

println("="^70)
println("CONUS GROUNDWATER STORAGE ANALYSIS")
println("="^70)
println("Package: GWLevel_USA v0.1.0")
println("Baseline: $(cfg.baseline_start)-$(cfg.baseline_end)")
println()

# Load harmonized data
data = load_harmonized_data(HARMONIZED_CSV; cfg=cfg)

# Filter by season
data = filter_by_season(data, cfg.valid_seasons)

println("Loaded $(nrow(data)) records from $(length(unique(data.WellID))) wells")
println("Aquifers: $(length(unique(data.Aquifer)))")

# =============================================================================
# Load Aquifer Areas
# =============================================================================

# Pre-computed areas (km²) from EPA Principal Aquifer boundaries
# Computed in EPSG:5070 (NAD83/Conus Albers Equal-Area)
aquifer_areas = Dict{String, Float64}(
    "N_Great_Plains" => 670433.0,
    "Great_lakes" => 653488.0,
    "High_Plains" => 505486.0,
    "Basin_and_Range_Basin_Fill" => 429865.0,
    "Piedmont_and_Blue_Ridge" => 376474.0,
    "U_Colorado" => 369858.0,
    "N_High_Plains" => 248320.0,
    "NE_Glaciers" => 244101.0,
    "Floridan" => 214839.0,
    "Arizona_Alluvials" => 209945.0,
    "C_S_High_Plains" => 203400.0,
    "Mississippi_Embayment" => 199098.0,
    "Snake" => 186489.0,
    "Pennsylvania" => 159559.0,
    "Cambrian_Ordovician" => 152207.0,
    "Coastal_Lowland" => 143031.0,
    "N_Atlantic_Coastal_Plain" => 134847.0,
    "Ozark_Plateau" => 134633.0,
    "Texas_Gulf_Coast" => 126522.0,
    "SE_Coastal_Plain" => 125424.0,
    "Rio_Grande" => 112813.0,
    "Columbia_Plateau_Basaltic_rock" => 110957.0,
    "Carrizo_Wilcox" => 104887.0,
    "Trinity" => 98266.0,
    "Edwards_Plateau" => 98126.0,
    "Mississippi_River_Valley_Alluvial" => 94635.0,
    "San_Joaquin" => 83272.0,
    "Sacramento" => 72283.0,
    "Columbia_Plateau_Basin_fill" => 62683.0,
    "Central_Valley" => 62421.0
)

# =============================================================================
# Run Analysis
# =============================================================================

println()
println("Running storage analysis...")
println()

results = run_storage_analysis(data; cfg=cfg, aquifer_areas=aquifer_areas)

# =============================================================================
# Print Summary
# =============================================================================

print_summary(results)

# =============================================================================
# Export Results
# =============================================================================

mkpath(OUTPUT_DIR)
export_results(results, OUTPUT_DIR)

println()
println("Done! Results saved to: $OUTPUT_DIR")

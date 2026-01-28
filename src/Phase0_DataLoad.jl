# =============================================================================
# Phase 0: Data Download and Loading
# =============================================================================
#
# Download and load raw groundwater data from multiple sources:
#   - NWIS (USGS) - via dataRetrieval R package
#   - California DWR - CSV files
#   - Texas TWDB - RDS file
#
# Download requires R with dataRetrieval package installed.
# Loading done in Julia for speed.
#
# Author: Dr. Ashraf Rateb, Bureau of Economic Geology, UT Austin
# =============================================================================

using CSV
using DataFrames
using Dates
using RData  # For reading .RData files

export download_nwis_data, download_all_data
export load_nwis_data, load_california_data, load_texas_data
export load_all_raw_data, harmonize_data

# CONUS states
const CONUS_STATES = [
    "AL", "AZ", "AR", "CA", "CO", "CT", "DE", "FL", "GA", "ID",
    "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD", "MA", "MI",
    "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY",
    "NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN",
    "TX", "UT", "VT", "VA", "WA", "WV", "WI", "WY"
]

"""
    download_nwis_data(output_dir::String; states=CONUS_STATES)

Download NWIS groundwater level data for specified states using R dataRetrieval.
Requires R with dataRetrieval and dplyr packages installed.

# Arguments
- `output_dir`: Directory to save .RData files
- `states`: Vector of state codes (default: all 48 CONUS states)

# Example
```julia
download_nwis_data("data/raw/NWIS")
download_nwis_data("data/raw/NWIS"; states=["TX", "CA", "AZ"])
```
"""
function download_nwis_data(output_dir::String; states=CONUS_STATES)
    mkpath(output_dir)

    # R script for downloading
    r_script = """
    library(dataRetrieval)
    library(dplyr)

    safe_bind <- function(df_list) {
      if (length(df_list) == 0) return(NULL)
      df_list <- lapply(df_list, function(df) { df[] <- lapply(df, as.character); df })
      bind_rows(df_list)
    }

    get_site_info_chunked <- function(site_nos, chunk_size = 300) {
      n <- length(site_nos)
      if (n == 0) return(NULL)
      all_info <- list()
      for (i in seq(1, n, by = chunk_size)) {
        end_i <- min(i + chunk_size - 1, n)
        chunk <- site_nos[i:end_i]
        info <- tryCatch(readNWISsite(siteNumbers = chunk), error = function(e) NULL)
        if (!is.null(info) && nrow(info) > 0) all_info[[length(all_info) + 1]] <- info
      }
      if (length(all_info) == 0) return(NULL)
      safe_bind(all_info)
    }

    get_gwl_chunked <- function(site_nos, chunk_size = 2000) {
      n <- length(site_nos)
      if (n == 0) return(NULL)
      all_data <- list()
      for (i in seq(1, n, by = chunk_size)) {
        end_i <- min(i + chunk_size - 1, n)
        chunk <- site_nos[i:end_i]
        data <- tryCatch(readNWISgwl(siteNumbers = chunk), error = function(e) {
          Sys.sleep(3); tryCatch(readNWISgwl(siteNumbers = chunk), error = function(e2) NULL)
        })
        if (!is.null(data) && nrow(data) > 0) all_data[[length(all_data) + 1]] <- data
        Sys.sleep(1)
      }
      if (length(all_data) == 0) return(NULL)
      safe_bind(all_data)
    }

    download_state <- function(state, output_dir) {
      tryCatch({
        sites <- whatNWISsites(stateCd = state, parameterCd = "72019")
        if (nrow(sites) == 0) return(list(state=state, n=0, status="no_sites"))

        n_sites <- nrow(sites)
        if (n_sites > 5000) {
          data <- get_gwl_chunked(sites\$site_no)
        } else {
          data <- readNWISgwl(siteNumbers = sites\$site_no)
        }

        if (is.null(data) || nrow(data) == 0) return(list(state=state, n=0, status="no_data"))

        unique_sites <- unique(data\$site_no)
        site_info <- get_site_info_chunked(unique_sites)

        if (!is.null(site_info)) {
          keep_cols <- c("site_no", "dec_lat_va", "dec_long_va", "well_depth_va",
                         "nat_aqfr_cd", "aqfr_cd", "aqfr_type_cd", "alt_va", "alt_datum_cd")
          keep_cols <- keep_cols[keep_cols %in% names(site_info)]
          site_info <- site_info[, keep_cols, drop = FALSE]
          data <- data %>% left_join(site_info, by = "site_no")
        }

        fname <- file.path(output_dir, paste0("NWIS_", state, "_", Sys.Date(), ".RData"))
        state_data <- data
        save(state_data, file = fname)
        return(list(state=state, n=nrow(data), status="ok"))
      }, error = function(e) {
        return(list(state=state, n=0, status=paste("error:", conditionMessage(e))))
      })
    }
    """

    println("=" ^ 60)
    println("DOWNLOADING NWIS GROUNDWATER DATA")
    println("=" ^ 60)
    println("States: $(length(states))")
    println("Output: $output_dir")
    println()

    for state in states
        print("[$state] ")

        # Write and run R script for this state
        state_script = r_script * """
        result <- download_state("$state", "$output_dir")
        cat(result\$state, result\$n, result\$status, sep="\\t")
        """

        try
            result = read(`Rscript -e $state_script`, String)
            parts = split(strip(result), '\t')
            if length(parts) >= 3
                n_records = tryparse(Int, parts[2])
                status = parts[3]
                if status == "ok" && n_records !== nothing
                    println("$n_records records")
                else
                    println(status)
                end
            else
                println("done")
            end
        catch e
            println("ERROR: $e")
        end
    end

    # Summary
    files = filter(f -> endswith(f, ".RData"), readdir(output_dir))
    println()
    println("=" ^ 60)
    println("Downloaded $(length(files)) state files to: $output_dir")
    println("=" ^ 60)

    return output_dir
end

"""
    download_all_data(base_dir::String)

Download all groundwater data sources (NWIS for all 48 CONUS states).

# Example
```julia
download_all_data("data/raw")
```
"""
function download_all_data(base_dir::String)
    nwis_dir = joinpath(base_dir, "NWIS")
    download_nwis_data(nwis_dir)
    return base_dir
end

"""
    load_nwis_data(data_dir::String) -> DataFrame

Load all NWIS .RData files from directory.
Returns standardized DataFrame with columns:
  WellID, Latitude, Longitude, MeasurementDate, DepthToWater_ft,
  WellDepth_ft, Altitude_ft, NatAquiferCode, LocalAquiferCode, Source
"""
function load_nwis_data(data_dir::String)
    println("Loading NWIS data from: $data_dir")

    all_data = DataFrame[]
    files = filter(f -> endswith(f, ".RData"), readdir(data_dir))

    for f in files
        try
            path = joinpath(data_dir, f)
            rdata = RData.load(path)

            # RData files contain 'state_data' variable
            if haskey(rdata, "state_data")
                df = DataFrame(rdata["state_data"])

                # Extract state from filename
                state = replace(f, r"NWIS_([A-Z]{2})_.*" => s"\1")

                # Standardize column names (NWIS format from dataRetrieval)
                rename_cols = Dict(
                    :site_no => :WellID,
                    :lev_dt => :MeasurementDate,
                    :lev_va => :DepthToWater_ft,
                    :dec_lat_va => :Latitude,
                    :dec_long_va => :Longitude,
                    :well_depth_va => :WellDepth_ft,
                    :alt_va => :Altitude_ft,
                    :nat_aqfr_cd => :NatAquiferCode,
                    :aqfr_cd => :LocalAquiferCode
                )

                for (old, new) in rename_cols
                    if hasproperty(df, old)
                        rename!(df, old => new)
                    end
                end

                df[!, :Source] .= "NWIS"
                df[!, :State] .= state
                push!(all_data, df)
                println("  $state: $(nrow(df)) records")
            end
        catch e
            @warn "Error loading $f: $e"
        end
    end

    if isempty(all_data)
        return DataFrame()
    end

    combined = vcat(all_data...; cols=:union)
    println("  TOTAL: $(nrow(combined)) records from $(length(files)) states")
    return combined
end

"""
    load_california_data(data_dir::String) -> DataFrame

Load California DWR groundwater data from CSV files.
"""
function load_california_data(data_dir::String)
    println("Loading California data from: $data_dir")

    # Load stations for coordinates
    stations_file = joinpath(data_dir, "gwl-stations.csv")
    if !isfile(stations_file)
        @warn "Stations file not found: $stations_file"
        return DataFrame()
    end

    stations = CSV.read(stations_file, DataFrame)

    # Load measurements (monthly is more complete)
    data_file = joinpath(data_dir, "gwl-monthly.csv")
    if !isfile(data_file)
        data_file = joinpath(data_dir, "gwl-daily.csv")
    end

    if !isfile(data_file)
        @warn "No measurement data found in $data_dir"
        return DataFrame()
    end

    df = CSV.read(data_file, DataFrame)

    # Merge with station coordinates
    if hasproperty(stations, :STATION) && hasproperty(df, :STATION)
        df = leftjoin(df,
                      select(stations, :STATION, :LATITUDE, :LONGITUDE);
                      on=:STATION)
    end

    # Standardize columns
    result = DataFrame(
        WellID = string.(df.STATION),
        Latitude = get(df, :LATITUDE, missing),
        Longitude = get(df, :LONGITUDE, missing),
        MeasurementDate = df.MSMT_DATE,
        DepthToWater_ft = df.GSE_WSE,  # Ground surface - water surface = depth
        Source = fill("CA_DWR", nrow(df))
    )

    println("  Loaded $(nrow(result)) records")
    return result
end

"""
    load_texas_data(data_dir::String) -> DataFrame

Load Texas TWDB groundwater data.
"""
function load_texas_data(data_dir::String)
    println("Loading Texas data from: $data_dir")

    # Look for RDS or RData file
    files = readdir(data_dir)
    rds_files = filter(f -> endswith(f, ".rds") || endswith(f, ".RData"), files)

    if isempty(rds_files)
        @warn "No Texas data files found in $data_dir"
        return DataFrame()
    end

    path = joinpath(data_dir, first(rds_files))

    try
        rdata = load(path)
        # Get the first dataframe from the RData
        df = DataFrame(first(values(rdata)))

        # Standardize columns (Texas format)
        rename_map = Dict(
            :StateWellNumber => :WellID,
            :LatitudeDD => :Latitude,
            :LongitudeDD => :Longitude,
            :DepthFromLSD => :DepthToWater_ft,
            :MeasurementDate => :MeasurementDate
        )

        for (old, new) in rename_map
            if hasproperty(df, old)
                rename!(df, old => new)
            end
        end

        df[!, :Source] .= "TX_TWDB"
        println("  Loaded $(nrow(df)) records")
        return df

    catch e
        @warn "Error loading Texas data: $e"
        return DataFrame()
    end
end

"""
    load_all_raw_data(base_dir::String) -> DataFrame

Load all raw data sources and combine.
"""
function load_all_raw_data(base_dir::String)
    println("=" ^ 60)
    println("LOADING RAW GROUNDWATER DATA")
    println("=" ^ 60)

    all_data = DataFrame[]

    # NWIS
    nwis_dir = joinpath(base_dir, "NWIS")
    if isdir(nwis_dir)
        nwis = load_nwis_data(nwis_dir)
        if nrow(nwis) > 0
            push!(all_data, nwis)
        end
    end

    # California
    ca_dir = joinpath(base_dir, "CA")
    if isdir(ca_dir)
        ca = load_california_data(ca_dir)
        if nrow(ca) > 0
            push!(all_data, ca)
        end
    end

    # Texas
    tx_dir = joinpath(base_dir, "TX")
    if isdir(tx_dir)
        tx = load_texas_data(tx_dir)
        if nrow(tx) > 0
            push!(all_data, tx)
        end
    end

    if isempty(all_data)
        error("No data loaded!")
    end

    combined = vcat(all_data...; cols=:union)
    println("\n" * "=" ^ 60)
    println("TOTAL: $(nrow(combined)) records from $(length(unique(combined.WellID))) wells")
    println("=" ^ 60)

    return combined
end

"""
    harmonize_data(df::DataFrame) -> DataFrame

Clean and harmonize the combined data:
- Convert units to meters
- Parse dates and add Year, Month, Season columns
- Filter to valid CONUS bounds
- Remove obvious errors
- Remove duplicates
"""
function harmonize_data(df::DataFrame)
    println("\nHarmonizing data...")
    n_start = nrow(df)

    result = copy(df)

    # Ensure WellID is string
    result.WellID = string.(result.WellID)

    # Parse MeasurementDate if it's a string
    if hasproperty(result, :MeasurementDate)
        if eltype(result.MeasurementDate) <: AbstractString
            result.MeasurementDate = Date.(result.MeasurementDate, dateformat"yyyy-mm-dd")
        end
    end

    # Convert depth to meters (handle string columns from safe_bind)
    if hasproperty(result, :DepthToWater_ft)
        depth_ft = result.DepthToWater_ft
        if eltype(depth_ft) <: AbstractString
            depth_ft = tryparse.(Float64, depth_ft)
        end
        result[!, :DepthToWater_m] = coalesce.(depth_ft, NaN) .* 0.3048
    end

    # Convert coordinates if strings
    if hasproperty(result, :Latitude) && eltype(result.Latitude) <: AbstractString
        result.Latitude = tryparse.(Float64, result.Latitude)
    end
    if hasproperty(result, :Longitude) && eltype(result.Longitude) <: AbstractString
        result.Longitude = tryparse.(Float64, result.Longitude)
    end

    # Add temporal columns
    if hasproperty(result, :MeasurementDate)
        result[!, :Year] = year.(result.MeasurementDate)
        result[!, :Month] = month.(result.MeasurementDate)

        # Add season
        result[!, :Season] = map(result.Month) do m
            if m in [12, 1, 2]
                "Winter"
            elseif m in [3, 4, 5]
                "Spring"
            elseif m in [6, 7, 8]
                "Summer"
            else
                "Fall"
            end
        end
    end

    # Filter CONUS bounds
    valid_lat = coalesce.(result.Latitude .>= 24, false) .&
                coalesce.(result.Latitude .<= 50, false)
    valid_lon = coalesce.(result.Longitude .>= -125, false) .&
                coalesce.(result.Longitude .<= -66, false)
    result = result[valid_lat .& valid_lon, :]
    println("  After CONUS filter: $(nrow(result)) records")

    # Filter valid depths (0 to 500m)
    if hasproperty(result, :DepthToWater_m)
        valid_depth = coalesce.(result.DepthToWater_m .> 0, false) .&
                      coalesce.(result.DepthToWater_m .< 500, false)
        result = result[valid_depth, :]
        println("  After depth filter: $(nrow(result)) records")
    end

    # Remove duplicates
    result = unique(result, [:WellID, :MeasurementDate])

    println("  Final: $(nrow(result)) records (removed $(n_start - nrow(result)) invalid/duplicate)")
    return result
end

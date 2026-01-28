#!/usr/bin/env Rscript
# =============================================================================
# Phase 0A: Download NWIS Groundwater Data
# =============================================================================
#
# Downloads groundwater level data from USGS National Water Information System
# for all CONUS states using the dataRetrieval package.
#
# Usage:
#   Rscript 00_download_nwis.R [output_dir]
#
# Output:
#   - NWIS_<STATE>_<DATE>.rds for each state
#   - NWIS_ALL_<DATE>.rds combined file
#
# Author: Dr. Ashraf Rateb, Bureau of Economic Geology, UT Austin
# =============================================================================

suppressPackageStartupMessages({
  library(dataRetrieval)
  library(dplyr)
  library(lubridate)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

CONUS_STATES <- c(
  "AL", "AZ", "AR", "CA", "CO", "CT", "DE", "FL", "GA", "ID",
  "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD", "MA", "MI",
  "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY",
  "NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN",
  "TX", "UT", "VT", "VA", "WA", "WV", "WI", "WY"
)

# =============================================================================
# FUNCTIONS
# =============================================================================

#' Download NWIS groundwater data for a single state
#'
#' @param state_abbr Two-letter state abbreviation
#' @param output_dir Directory to save output
#' @return Data frame with groundwater levels, or NULL on error
download_nwis_state <- function(state_abbr, output_dir = ".") {
  cat(sprintf("  [%s] Downloading...\n", state_abbr))

  tryCatch({
    # Get all groundwater monitoring sites
    sites <- whatNWISsites(
      stateCd = state_abbr,
      parameterCd = "72019"  # Depth to water level, feet below land surface
    )

    if (is.null(sites) || nrow(sites) == 0) {
      cat(sprintf("  [%s] No sites found\n", state_abbr))
      return(NULL)
    }

    cat(sprintf("  [%s] Found %d sites, fetching data...\n", state_abbr, nrow(sites)))

    # Download groundwater levels
    gw_data <- readNWISgwl(siteNumbers = sites$site_no)

    if (is.null(gw_data) || nrow(gw_data) == 0) {
      cat(sprintf("  [%s] No data retrieved\n", state_abbr))
      return(NULL)
    }

    # Merge with site metadata
    result <- gw_data %>%
      left_join(
        sites %>% select(site_no, dec_lat_va, dec_long_va, station_nm,
                         nat_aqfr_cd, aqfr_cd, well_depth_va),
        by = "site_no"
      ) %>%
      transmute(
        WellID = site_no,
        WellName = station_nm,
        Latitude = dec_lat_va,
        Longitude = dec_long_va,
        MeasurementDate = as.Date(lev_dt),
        DepthToWater_ft = lev_va,
        WellDepth_ft = well_depth_va,
        Aquifer_National = nat_aqfr_cd,
        Aquifer_Local = aqfr_cd,
        State = state_abbr,
        Source = "NWIS"
      )

    cat(sprintf("  [%s] Downloaded %d records from %d wells\n",
                state_abbr, nrow(result), length(unique(result$WellID))))

    # Save state file
    out_file <- file.path(output_dir, sprintf("NWIS_%s_%s.rds", state_abbr, Sys.Date()))
    saveRDS(result, out_file)
    cat(sprintf("  [%s] Saved: %s\n", state_abbr, basename(out_file)))

    return(result)

  }, error = function(e) {
    cat(sprintf("  [%s] ERROR: %s\n", state_abbr, conditionMessage(e)))
    return(NULL)
  })
}

#' Download NWIS groundwater data for all CONUS states
#'
#' @param output_dir Directory to save output files
#' @param states Vector of state abbreviations (default: all CONUS)
#' @return Combined data frame
download_nwis_all <- function(output_dir = ".", states = CONUS_STATES) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cat("=========================================================\n")
  cat("NWIS GROUNDWATER DATA DOWNLOAD\n")
  cat("=========================================================\n")
  cat(sprintf("Date: %s\n", Sys.Date()))
  cat(sprintf("States: %d\n", length(states)))
  cat(sprintf("Output: %s\n\n", output_dir))

  # Download each state
  all_data <- list()
  for (state in states) {
    result <- download_nwis_state(state, output_dir)
    if (!is.null(result)) {
      all_data[[state]] <- result
    }
  }

  # Combine all states
  if (length(all_data) > 0) {
    combined <- bind_rows(all_data)

    # Save combined file
    out_file <- file.path(output_dir, sprintf("NWIS_ALL_%s.rds", Sys.Date()))
    saveRDS(combined, out_file)

    cat("\n=========================================================\n")
    cat("DOWNLOAD COMPLETE\n")
    cat("=========================================================\n")
    cat(sprintf("Total records: %s\n", format(nrow(combined), big.mark = ",")))
    cat(sprintf("Total wells: %s\n", format(length(unique(combined$WellID)), big.mark = ",")))
    cat(sprintf("States: %d\n", length(all_data)))
    cat(sprintf("Combined file: %s\n", out_file))

    return(combined)
  } else {
    cat("ERROR: No data downloaded\n")
    return(NULL)
  }
}

# =============================================================================
# MAIN
# =============================================================================

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  output_dir <- if (length(args) >= 1) args[1] else "data/raw/NWIS"

  download_nwis_all(output_dir)
}

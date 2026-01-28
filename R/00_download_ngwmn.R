#!/usr/bin/env Rscript
# =============================================================================
# Phase 0B: Download NGWMN Groundwater Data
# =============================================================================
#
# Downloads groundwater level data from the National Groundwater Monitoring
# Network (NGWMN) using the dataRetrieval package.
#
# Usage:
#   Rscript 00_download_ngwmn.R [output_dir]
#
# Output:
#   - NGWMN_sites_<DATE>.rds (site metadata)
#   - NGWMN_ALL_<DATE>.rds (observations)
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

# CONUS bounding box (minLat, minLon, maxLat, maxLon)
CONUS_BBOX <- c(24, -125, 50, -66)

# =============================================================================
# FUNCTIONS
# =============================================================================

#' Download NGWMN site metadata
#'
#' @param bbox Bounding box (minLat, minLon, maxLat, maxLon)
#' @return Data frame with site information
download_ngwmn_sites <- function(bbox = CONUS_BBOX) {
  cat("Fetching NGWMN site list...\n")

  tryCatch({
    sites <- readNGWMNdata(
      service = "featureOfInterest",
      bbox = bbox
    )

    if (is.null(sites) || nrow(sites) == 0) {
      cat("No NGWMN sites found in bounding box\n")
      return(NULL)
    }

    cat(sprintf("Found %d NGWMN sites\n", nrow(sites)))
    return(sites)

  }, error = function(e) {
    cat(sprintf("ERROR fetching sites: %s\n", conditionMessage(e)))
    return(NULL)
  })
}

#' Download observations for a single NGWMN site
#'
#' @param site_id Site identifier
#' @param lat Site latitude
#' @param lon Site longitude
#' @return Data frame with observations, or NULL on error
download_ngwmn_site_obs <- function(site_id, lat, lon) {
  tryCatch({
    data <- readNGWMNdata(
      service = "observation",
      siteNumbers = site_id
    )

    if (is.null(data) || nrow(data) == 0) {
      return(NULL)
    }

    data$Latitude <- lat
    data$Longitude <- lon
    return(data)

  }, error = function(e) {
    return(NULL)
  })
}

#' Download all NGWMN groundwater data for CONUS
#'
#' @param output_dir Directory to save output files
#' @param bbox Bounding box for site query
#' @return Combined data frame
download_ngwmn_all <- function(output_dir = ".", bbox = CONUS_BBOX) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cat("=========================================================\n")
  cat("NGWMN GROUNDWATER DATA DOWNLOAD\n")
  cat("=========================================================\n")
  cat(sprintf("Date: %s\n", Sys.Date()))
  cat(sprintf("Bounding box: %.1f, %.1f, %.1f, %.1f\n",
              bbox[1], bbox[2], bbox[3], bbox[4]))
  cat(sprintf("Output: %s\n\n", output_dir))

  # Get site list
  sites <- download_ngwmn_sites(bbox)
  if (is.null(sites)) return(NULL)

  # Save site metadata
  sites_file <- file.path(output_dir, sprintf("NGWMN_sites_%s.rds", Sys.Date()))
  saveRDS(sites, sites_file)
  cat(sprintf("Saved site metadata: %s\n\n", basename(sites_file)))

  # Download observations for each site
  cat("Downloading observations...\n")
  n_sites <- nrow(sites)
  all_data <- list()
  success_count <- 0

  for (i in seq_len(n_sites)) {
    if (i %% 100 == 0 || i == n_sites) {
      cat(sprintf("  Progress: %d/%d sites (%.1f%%)\n",
                  i, n_sites, 100 * i / n_sites))
    }

    result <- download_ngwmn_site_obs(
      site_id = sites$site[i],
      lat = sites$dec_lat_va[i],
      lon = sites$dec_lon_va[i]
    )

    if (!is.null(result) && nrow(result) > 0) {
      all_data[[length(all_data) + 1]] <- result
      success_count <- success_count + 1
    }
  }

  if (length(all_data) == 0) {
    cat("ERROR: No observation data retrieved\n")
    return(NULL)
  }

  # Combine and standardize
  combined <- bind_rows(all_data) %>%
    transmute(
      WellID = site,
      Latitude = Latitude,
      Longitude = Longitude,
      MeasurementDate = as.Date(date),
      DepthToWater_ft = value,  # Note: Check units from NGWMN
      Source = "NGWMN"
    )

  # Save combined file
  out_file <- file.path(output_dir, sprintf("NGWMN_ALL_%s.rds", Sys.Date()))
  saveRDS(combined, out_file)

  cat("\n=========================================================\n")
  cat("DOWNLOAD COMPLETE\n")
  cat("=========================================================\n")
  cat(sprintf("Sites queried: %d\n", n_sites))
  cat(sprintf("Sites with data: %d\n", success_count))
  cat(sprintf("Total records: %s\n", format(nrow(combined), big.mark = ",")))
  cat(sprintf("Total wells: %s\n", format(length(unique(combined$WellID)), big.mark = ",")))
  cat(sprintf("Output file: %s\n", out_file))

  return(combined)
}

# =============================================================================
# MAIN
# =============================================================================

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  output_dir <- if (length(args) >= 1) args[1] else "data/raw/NGWMN"

  download_ngwmn_all(output_dir)
}

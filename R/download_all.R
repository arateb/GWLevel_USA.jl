#!/usr/bin/env Rscript
# =============================================================================
# Phase 0: Download All Groundwater Data
# =============================================================================
#
# Download fresh data from NWIS with aquifer info using dataRetrieval.
# See: https://waterdata.usgs.gov/blog/dataretrieval/
#
# Usage: Rscript download_all.R
#
# Author: Dr. Ashraf Rateb
# =============================================================================

library(dataRetrieval)
library(dplyr)

OUTPUT_DIR <- "data/raw/NWIS"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

CONUS_STATES <- c(
  "AL", "AZ", "AR", "CA", "CO", "CT", "DE", "FL", "GA", "ID",
  "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD", "MA", "MI",
  "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY",
  "NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN",
  "TX", "UT", "VT", "VA", "WA", "WV", "WI", "WY"
)

# Helper: get site info in chunks to avoid HTTP 400
get_site_info_chunked <- function(site_nos, chunk_size = 300) {
  n <- length(site_nos)
  if (n == 0) return(NULL)

  all_info <- list()
  for (i in seq(1, n, by = chunk_size)) {
    end_i <- min(i + chunk_size - 1, n)
    chunk <- site_nos[i:end_i]
    info <- tryCatch(
      readNWISsite(siteNumbers = chunk),
      error = function(e) NULL
    )
    if (!is.null(info) && nrow(info) > 0) {
      all_info[[length(all_info) + 1]] <- info
    }
  }

  if (length(all_info) == 0) return(NULL)
  bind_rows(all_info)
}

cat("Downloading NWIS groundwater data for", length(CONUS_STATES), "states\n")
cat("Output:", OUTPUT_DIR, "\n\n")

for (state in CONUS_STATES) {
  cat(sprintf("[%s] ", state))

  tryCatch({
    # Get sites with GW level parameter 72019 (depth to water)
    sites <- whatNWISsites(stateCd = state, parameterCd = "72019")
    if (nrow(sites) == 0) { cat("no sites\n"); next }

    cat(sprintf("%d sites, ", nrow(sites)))

    # Get groundwater levels
    data <- readNWISgwl(siteNumbers = sites$site_no)
    if (nrow(data) == 0) { cat("no data\n"); next }

    # Get full site info in chunks
    unique_sites <- unique(data$site_no)
    site_info <- get_site_info_chunked(unique_sites)

    if (!is.null(site_info)) {
      # Select columns that exist
      keep_cols <- c("site_no", "dec_lat_va", "dec_long_va",
                     "well_depth_va", "nat_aqfr_cd", "aqfr_cd",
                     "aqfr_type_cd", "alt_va", "alt_datum_cd")
      keep_cols <- keep_cols[keep_cols %in% names(site_info)]
      site_info <- site_info[, keep_cols, drop = FALSE]

      # Merge GW levels with site info
      data <- data %>%
        left_join(site_info, by = "site_no")
    }

    # Save
    fname <- file.path(OUTPUT_DIR, paste0("NWIS_", state, "_", Sys.Date(), ".RData"))
    state_data <- data
    save(state_data, file = fname)
    cat(nrow(data), "records\n")

  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
  })
}

cat("\nDone! Files saved to:", OUTPUT_DIR, "\n")

# Summary
files <- list.files(OUTPUT_DIR, pattern = "*.RData", full.names = TRUE)
cat("\nSummary:\n")
cat("  Total files:", length(files), "\n")
total_records <- 0
for (f in files) {
  load(f)
  total_records <- total_records + nrow(state_data)
}
cat("  Total records:", formatC(total_records, big.mark = ","), "\n")

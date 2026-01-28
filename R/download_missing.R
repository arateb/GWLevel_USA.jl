#!/usr/bin/env Rscript
# Download missing large states with chunked download

library(dataRetrieval)
library(dplyr)

OUTPUT_DIR <- "data/raw/NWIS"

# Chunk helper for site info
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
  bind_rows(all_info)
}

# Chunk helper for GW levels (for mega-states)
get_gwl_chunked <- function(site_nos, chunk_size = 2000) {
  n <- length(site_nos)
  if (n == 0) return(NULL)
  all_data <- list()
  for (i in seq(1, n, by = chunk_size)) {
    end_i <- min(i + chunk_size - 1, n)
    chunk <- site_nos[i:end_i]
    cat(sprintf("  batch %d-%d of %d sites... ", i, end_i, n))
    data <- tryCatch(readNWISgwl(siteNumbers = chunk), error = function(e) {
      cat("error, retrying in 5s... ")
      Sys.sleep(5)
      tryCatch(readNWISgwl(siteNumbers = chunk), error = function(e2) NULL)
    })
    if (!is.null(data) && nrow(data) > 0) {
      all_data[[length(all_data) + 1]] <- data
      cat(nrow(data), "records\n")
    } else {
      cat("no data\n")
    }
    Sys.sleep(1)  # Be nice to the server
  }
  if (length(all_data) == 0) return(NULL)
  bind_rows(all_data)
}

MISSING <- c("AZ", "CA", "KS", "KY", "MN", "MS", "PA", "WA")

for (state in MISSING) {
  cat(sprintf("\n[%s] ", state))
  tryCatch({
    sites <- whatNWISsites(stateCd = state, parameterCd = "72019")
    if (nrow(sites) == 0) { cat("no sites\n"); next }
    n_sites <- nrow(sites)
    cat(sprintf("%d sites\n", n_sites))

    # For mega-states, chunk the GW level download
    if (n_sites > 5000) {
      cat("  Large state - downloading in chunks...\n")
      data <- get_gwl_chunked(sites$site_no)
    } else {
      data <- readNWISgwl(siteNumbers = sites$site_no)
    }

    if (is.null(data) || nrow(data) == 0) { cat("  no data\n"); next }
    cat(sprintf("  Total: %d records\n", nrow(data)))

    # Get site info in chunks
    unique_sites <- unique(data$site_no)
    cat(sprintf("  Getting site info for %d unique wells...\n", length(unique_sites)))
    site_info <- get_site_info_chunked(unique_sites)

    if (!is.null(site_info)) {
      keep_cols <- c("site_no", "dec_lat_va", "dec_long_va", "well_depth_va",
                     "nat_aqfr_cd", "aqfr_cd", "aqfr_type_cd", "alt_va", "alt_datum_cd")
      keep_cols <- keep_cols[keep_cols %in% names(site_info)]
      site_info <- site_info[, keep_cols, drop = FALSE]
      data <- data %>% left_join(site_info, by = "site_no")
    }

    fname <- file.path(OUTPUT_DIR, paste0("NWIS_", state, "_", Sys.Date(), ".RData"))
    state_data <- data
    save(state_data, file = fname)
    cat(sprintf("  SAVED: %s (%d records)\n", fname, nrow(data)))
  }, error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
}

cat("\n\nDone!\n")
cat("Total files:", length(list.files(OUTPUT_DIR, pattern = "*.RData")), "\n")

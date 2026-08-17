################################################################################
#                                                                              #
#              SYNGAP1 Variant Viewer — Cache Builder                          #
#                                                                              #
#  Purpose: Standalone script that pre-builds both disk caches used by the     #
#           SVV Shiny app, independently of the Shiny runtime.                 #
#                                                                              #
#  Caches built:                                                               #
#    1. ClinVar data  — cache/clinvar/                                         #
#         clinvar_data.rds     : parsed variant data frame                     #
#         last_updated.txt     : ISO timestamp of last successful fetch        #
#         Source: NCBI E-utilities (esearch + esummary)                        #
#         Expiry: 7 days                                                       #
#                                                                              #
#    2. Transcript exon structure  — cache/transcript_structure/               #
#         ENST00000418600.rds            : exon data frame                     #
#         ENST00000418600.rds.timestamp  : ISO timestamp of last fetch         #
#         Source: Ensembl REST API                                             #
#         Expiry: 6 months (180 days)                                          #
#                                                                              #
#  Compatible with: app.R v3.2.0+                                          #
#  Note: As of v3.2.0 these functions are maintained solely here; they are no  #
#        longer duplicated in app.R.                                        #
#                                                                              #
#  Usage:                                                                      #
#    Rscript svv_build_cache.R          # from terminal                        #
#    source("svv_build_cache.R")        # from R console or RStudio            #
#                                                                              #
#  Diagnostic intent:                                                          #
#    Running this script on the server (or locally) confirms whether:          #
#    (a) outbound HTTP to NCBI and Ensembl is permitted, and                   #
#    (b) the process has write permission to create the cache directories.     #
#    All status messages are printed to the console with [OK] / [WARN] /       #
#    [FAIL] prefixes so failures are easy to spot in a log or terminal.        #
#                                                                              #
################################################################################

################################################################################
# SECTION 1: PACKAGE DEPENDENCIES
################################################################################

library(httr)
library(jsonlite)

################################################################################
# SECTION 2: CLINVAR CONSTANTS
# ClinVar cache constants (also read by app.R via CLINVAR_CACHE_RDS)
################################################################################

CLINVAR_CACHE_RDS       <- "cache/clinvar/clinvar_data.rds"
CLINVAR_CACHE_TIMESTAMP <- "cache/clinvar/last_updated.txt"
CLINVAR_CACHE_DAYS      <- 7          # Re-fetch after this many days
CLINVAR_BATCH_SIZE      <- 200        # IDs per esummary call (NCBI recommended max)
CLINVAR_RATE_DELAY      <- 0.4        # Seconds between requests (stays under 3/sec)
NCBI_BASE               <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

################################################################################
# SECTION 3: TRANSCRIPT CONSTANTS
# Transcript cache constants (also read by app.R via TRANSCRIPT_CACHE_FILE)
################################################################################

TRANSCRIPT_ID          <- "ENST00000418600"
TRANSCRIPT_CACHE_DIR   <- "cache/transcript_structure"
TRANSCRIPT_CACHE_FILE  <- file.path(TRANSCRIPT_CACHE_DIR, paste0(TRANSCRIPT_ID, ".rds"))
TRANSCRIPT_CACHE_TS    <- paste0(TRANSCRIPT_CACHE_FILE, ".timestamp")
TRANSCRIPT_CACHE_DAYS  <- 180

################################################################################
# SECTION 4: CLINVAR HELPER FUNCTIONS
# Source: ClinVar cache-building functions (canonical home: svv_build_cache.R)
################################################################################

#' Check Whether the ClinVar Cache is Still Fresh
#'
#' Returns TRUE if a valid cache file exists and was written within
#' CLINVAR_CACHE_DAYS days. Returns FALSE otherwise.
clinvar_cache_is_fresh <- function() {
  if (!file.exists(CLINVAR_CACHE_RDS) ||
      !file.exists(CLINVAR_CACHE_TIMESTAMP)) return(FALSE)

  timestamp_str <- tryCatch(
    readLines(CLINVAR_CACHE_TIMESTAMP, n = 1, warn = FALSE),
    error = function(e) return("")
  )
  if (nchar(trimws(timestamp_str)) == 0) return(FALSE)

  last_updated <- tryCatch(
    as.POSIXct(timestamp_str, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC"),
    error = function(e) return(NA)
  )
  if (is.na(last_updated)) return(FALSE)

  age_days <- as.numeric(difftime(Sys.time(), last_updated, units = "days"))
  return(age_days < CLINVAR_CACHE_DAYS)
}

#' Fetch All SYNGAP1 ClinVar Variation IDs via esearch
#'
#' Queries NCBI esearch for all ClinVar entries associated with SYNGAP1.
#' Returns a character vector of variation IDs, or NULL on failure.
fetch_clinvar_ids <- function() {
  message("ClinVar: querying NCBI esearch for SYNGAP1 variation IDs...")

  url <- paste0(
    NCBI_BASE, "esearch.fcgi",
    "?db=clinvar",
    "&term=SYNGAP1%5Bgene%5D",
    "&retmax=10000",
    "&retmode=json"
  )

  message("  URL: ", url)

  resp <- tryCatch(
    GET(url, httr::timeout(30)),
    error = function(e) {
      warning("ClinVar esearch network error: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(resp)) return(NULL)

  if (http_error(resp)) {
    warning("ClinVar esearch HTTP error: status ", status_code(resp),
            " — ", content(resp, "text", encoding = "UTF-8"))
    return(NULL)
  }

  parsed <- tryCatch(
    fromJSON(content(resp, "text", encoding = "UTF-8")),
    error = function(e) {
      warning("ClinVar esearch JSON parse error: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(parsed)) return(NULL)

  ids <- parsed$esearchresult$idlist
  message("ClinVar: found ", length(ids), " variation IDs.")
  return(ids)
}

#' Parse a Single esummary Record into a One-Row Data Frame
#'
#' Extracts the fields we need from one ClinVar esummary JSON result object.
#' Any field that is absent or malformed in the JSON is set to NA gracefully.
#'
#' @param rec One list element from the esummary "result" object
#' @return Single-row data frame with all ClinVar metadata columns
parse_esummary_record <- function(rec) {

  safe <- function(x, default = NA_character_) {
    if (is.null(x) || length(x) == 0) return(default)
    v <- trimws(as.character(x[[1]]))
    if (length(v) == 0 || is.na(v) || v == "") return(default)
    v
  }

  vset         <- rec$variation_set[[1]]
  name         <- safe(vset$variation_name)
  gene         <- safe(rec$genes[[1]]$symbol)
  protein_chg  <- safe(rec$protein_change)
  accession    <- safe(rec$accession)
  variation_id <- safe(rec$uid)
  variant_type <- safe(rec$obj_type)
  canonical_spdi <- safe(vset$canonical_spdi)

  mol_conseq <- if (!is.null(rec$molecular_consequence_list) &&
                    length(rec$molecular_consequence_list) > 0) {
    safe(rec$molecular_consequence_list[[1]])
  } else { NA_character_ }

  allele_id <- NA_character_
  dbsnp_id  <- NA_character_

  germ      <- rec$germline_classification
  germ_sig  <- safe(germ$description)
  germ_date <- safe(germ$last_evaluated)
  germ_rev  <- safe(germ$review_status)

  condition <- if (!is.null(germ$trait_set) && length(germ$trait_set) > 0) {
    trait_names <- sapply(germ$trait_set, function(t) {
      if (is.character(t)) safe(t)
      else safe(t$trait_name)
    })
    trait_names <- trait_names[!is.na(trait_names)]
    if (length(trait_names) > 0) paste(trait_names, collapse = "; ") else NA_character_
  } else { NA_character_ }

  som       <- rec$clinical_impact_classification
  som_sig   <- safe(som$description)
  som_date  <- safe(som$last_evaluated)
  som_rev   <- safe(som$review_status)

  onco      <- rec$oncogenicity_classification
  onco_sig  <- safe(onco$description)
  onco_date <- safe(onco$last_evaluated)
  onco_rev  <- safe(onco$review_status)

  chr38   <- NA_character_
  start38 <- NA_real_
  end38   <- NA_real_

  locs <- vset$variation_loc
  if (!is.null(locs) && length(locs) > 0) {
    for (loc in locs) {
      if (!is.null(loc$assembly_name) && loc$assembly_name == "GRCh38") {
        chr38   <- safe(loc$chr)
        start38 <- suppressWarnings(as.numeric(safe(loc$start)))
        end38   <- suppressWarnings(as.numeric(safe(loc$stop)))
        # Note: start == stop for single-base variants; no +1 offset applied.
        # (Bug fixed v3.0.1 — see app.R parse_esummary_record comments)
        break
      }
    }
  }

  data.frame(
    Name                                          = name,
    Gene.s.                                       = gene,
    Protein.change                                = protein_chg,
    Condition.s.                                  = condition,
    Accession                                     = accession,
    VariationID                                   = variation_id,
    AlleleID.s.                                   = allele_id,
    dbSNP.ID                                      = dbsnp_id,
    Canonical.SPDI                                = canonical_spdi,
    Variant.type                                  = variant_type,
    Molecular.consequence                         = mol_conseq,
    Germline.classification                       = germ_sig,
    Germline.date.last.evaluated                  = germ_date,
    Germline.review.status                        = germ_rev,
    Somatic.clinical.impact                       = som_sig,
    Somatic.clinical.impact.date.last.evaluated   = som_date,
    Somatic.clinical.impact.review.status         = som_rev,
    Oncogenicity.classification                   = onco_sig,
    Oncogenicity.date.last.evaluated              = onco_date,
    Oncogenicity.review.status                    = onco_rev,
    clinvar_chr                                   = if (!is.na(chr38)) paste0("chr", chr38) else NA_character_,
    clinvar_start                                 = start38,
    clinvar_end                                   = end38,
    stringsAsFactors = FALSE
  )
}

#' Fetch Full ClinVar Metadata for a Vector of Variation IDs via esummary
#'
#' Calls NCBI esummary in batches of CLINVAR_BATCH_SIZE IDs, parses each
#' result record, and returns a single combined data frame.
#'
#' @param ids Character vector of ClinVar variation IDs from fetch_clinvar_ids()
#' @return Data frame of all variants, or NULL on complete failure
fetch_clinvar_summaries <- function(ids) {
  batches    <- split(ids, ceiling(seq_along(ids) / CLINVAR_BATCH_SIZE))
  n_batches  <- length(batches)
  all_rows   <- vector("list", n_batches)

  message("ClinVar: fetching summaries in ", n_batches,
          " batches of up to ", CLINVAR_BATCH_SIZE, " IDs each...")

  for (i in seq_along(batches)) {
    message("  Batch ", i, " / ", n_batches, "...")

    id_string <- paste(batches[[i]], collapse = ",")
    url <- paste0(
      NCBI_BASE, "esummary.fcgi",
      "?db=clinvar",
      "&id=", id_string,
      "&retmode=json"
    )

    resp <- tryCatch(GET(url), error = function(e) NULL)
    if (is.null(resp) || http_error(resp)) {
      warning("  esummary batch ", i, " failed — skipping.")
      Sys.sleep(CLINVAR_RATE_DELAY)
      next
    }

    parsed <- tryCatch(
      fromJSON(content(resp, "text", encoding = "UTF-8"), simplifyVector = FALSE),
      error = function(e) NULL
    )
    if (is.null(parsed) || is.null(parsed$result)) {
      warning("  esummary batch ", i, " returned unparseable JSON — skipping.")
      Sys.sleep(CLINVAR_RATE_DELAY)
      next
    }

    uids    <- parsed$result$uids
    records <- parsed$result[names(parsed$result) != "uids"]

    batch_rows <- lapply(records, function(rec) {
      tryCatch(parse_esummary_record(rec), error = function(e) NULL)
    })
    batch_rows <- Filter(Negate(is.null), batch_rows)

    if (length(batch_rows) > 0) {
      all_rows[[i]] <- do.call(rbind, batch_rows)
    }

    Sys.sleep(CLINVAR_RATE_DELAY)
  }

  combined <- do.call(rbind, Filter(Negate(is.null), all_rows))
  if (is.null(combined) || nrow(combined) == 0) return(NULL)

  message("ClinVar: successfully parsed ", nrow(combined), " variant records.")
  return(combined)
}

#' Load ClinVar Data — from Cache if Fresh, Otherwise Fetch from NCBI
#'
#' Decision logic:
#'   1. Cache fresh (< 7 days)  → load .rds instantly
#'   2. Cache absent or expired → fetch from NCBI, save to cache
#'   3. Fetch fails + stale cache exists → use stale cache with warning
#'   4. Fetch fails + no cache → return empty data frame
#'
#' @return Data frame ready for use by create_clinvar_gff3_data()
load_clinvar_data <- function() {

  if (clinvar_cache_is_fresh()) {
    message("ClinVar: loading from cache (< ", CLINVAR_CACHE_DAYS, " days old).")
    return(readRDS(CLINVAR_CACHE_RDS))
  }

  message("ClinVar: cache absent or expired — fetching from NCBI E-utilities...")

  ids  <- fetch_clinvar_ids()
  data <- if (!is.null(ids) && length(ids) > 0) fetch_clinvar_summaries(ids) else NULL

  if (!is.null(data) && nrow(data) > 0) {
    saveRDS(data, CLINVAR_CACHE_RDS)
    writeLines(
      format(Sys.time(), format = "%Y-%m-%dT%H:%M:%S", tz = "UTC"),
      CLINVAR_CACHE_TIMESTAMP
    )
    message("ClinVar: cache updated successfully.")
    return(data)
  }

  if (file.exists(CLINVAR_CACHE_RDS)) {
    warning(
      "ClinVar: NCBI fetch failed. Falling back to stale cache. ",
      "Data may be outdated."
    )
    return(readRDS(CLINVAR_CACHE_RDS))
  }

  warning(
    "ClinVar: NCBI fetch failed and no cache found. ",
    "ClinVar track will be empty."
  )
  return(data.frame())
}

################################################################################
# SECTION 5: TRANSCRIPT HELPER FUNCTIONS
# Source: Transcript cache-building functions (canonical home: svv_build_cache.R)
################################################################################

#' Fetch Transcript Exon Structure from Ensembl
#'
#' Retrieves the complete exon structure for a transcript via a single Ensembl
#' REST API call. Returns a tidy data frame of exons sorted in transcript order.
#'
#' @param transcript_id Ensembl transcript ID (default: ENST00000418600)
#' @return Data frame with columns: strand, exon_start, exon_end, exon_length
#'         Returns NULL on network or parse failure.
fetch_transcript_exons <- function(transcript_id = TRANSCRIPT_ID) {
  message("Ensembl: fetching exon structure for ", transcript_id, "...")

  server <- "https://rest.ensembl.org"
  ext    <- paste0("/lookup/id/", transcript_id, "?expand=1")
  url    <- paste0(server, ext)

  message("  URL: ", url)

  response <- GET(url, content_type("application/json"))

  if (http_error(response)) {
    message("Error fetching transcript structure. Status: ", status_code(response))
    return(NULL)
  }

  result <- fromJSON(content(response, as = "text", encoding = "UTF-8"))

  if (!"Exon" %in% names(result) || length(result$Exon) == 0) {
    message("No exon data in Ensembl response for ", transcript_id)
    return(NULL)
  }

  exons  <- result$Exon
  strand <- result$strand   # +1 or -1

  exon_df <- data.frame(
    strand      = strand,
    exon_start  = exons$start,
    exon_end    = exons$end,
    exon_length = exons$end - exons$start + 1L,
    stringsAsFactors = FALSE
  )

  exon_df <- exon_df[order(exon_df$exon_start, decreasing = (strand == -1L)), ]
  rownames(exon_df) <- NULL

  message("Ensembl: retrieved ", nrow(exon_df), " exons (strand ", strand, ")")
  return(exon_df)
}

#' Write Transcript Cache to Disk
#'
#' Saves the exon data frame and records the current UTC time as the fetch
#' timestamp. Both files must exist for the cache to be considered fresh.
write_transcript_cache <- function(exon_df) {
  saveRDS(exon_df, TRANSCRIPT_CACHE_FILE)
  writeLines(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S", tz = "UTC"),
    TRANSCRIPT_CACHE_TS
  )
  message("Transcript exon structure cached to ", TRANSCRIPT_CACHE_FILE,
          " (refreshes in ", TRANSCRIPT_CACHE_DAYS, " days)")
}

#' Check Whether the Transcript Cache is Still Fresh
#'
#' Returns TRUE if both the .rds and .timestamp files exist and the cache is
#' younger than TRANSCRIPT_CACHE_DAYS. Returns FALSE in any other case.
transcript_cache_is_fresh <- function() {
  if (!file.exists(TRANSCRIPT_CACHE_FILE) ||
      !file.exists(TRANSCRIPT_CACHE_TS)) return(FALSE)

  timestamp_str <- tryCatch(
    readLines(TRANSCRIPT_CACHE_TS, n = 1, warn = FALSE),
    error = function(e) ""
  )
  if (nchar(trimws(timestamp_str)) == 0) return(FALSE)

  last_updated <- tryCatch(
    as.POSIXct(timestamp_str, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC"),
    error = function(e) NA
  )
  if (is.na(last_updated)) return(FALSE)

  age_days <- as.numeric(difftime(Sys.time(), last_updated, units = "days"))
  return(age_days < TRANSCRIPT_CACHE_DAYS)
}

################################################################################
# SECTION 6: MAIN — BUILD BOTH CACHES WITH DIAGNOSTIC REPORTING
#
# This section mirrors the startup execution blocks from app.R but wraps
# each step in explicit [OK] / [WARN] / [FAIL] reporting so the output is
# readable at a glance, whether run interactively or piped to a log file.
################################################################################

cat("\n")
cat("================================================================================\n")
cat("  SVV Cache Builder\n")
cat("  Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("================================================================================\n\n")

# ── Step 1: Create cache directories ─────────────────────────────────────────

cat("--- Step 1: Creating cache directories ---\n")

dir_clinvar     <- "cache/clinvar"
dir_transcript  <- TRANSCRIPT_CACHE_DIR

for (d in c(dir_clinvar, dir_transcript)) {
  result <- tryCatch({
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    TRUE
  }, error = function(e) {
    cat("[FAIL] Could not create directory:", d, "\n")
    cat("       Reason:", conditionMessage(e), "\n")
    FALSE
  })

  if (result) {
    if (dir.exists(d)) {
      cat("[OK]   Directory ready:", d, "\n")
    } else {
      # dir.create() returned without error but dir still doesn't exist —
      # this is the write-permission failure mode on read-only deployments.
      cat("[FAIL] Directory does not exist after creation attempt:", d, "\n")
      cat("       This likely means the process lacks write permission here.\n")
    }
  }
}

cat("\n")

# ── Step 2: ClinVar cache ─────────────────────────────────────────────────────

cat("--- Step 2: ClinVar cache ---\n")

clinvar_data <- tryCatch({
  load_clinvar_data()
}, error = function(e) {
  cat("[FAIL] load_clinvar_data() threw an unexpected error:\n")
  cat("       ", conditionMessage(e), "\n")
  data.frame()
})

if (!is.null(clinvar_data) && nrow(clinvar_data) > 0) {
  cat("[OK]   ClinVar data loaded:", nrow(clinvar_data), "variants\n")
  if (file.exists(CLINVAR_CACHE_RDS)) {
    cat("[OK]   Cache file present:", CLINVAR_CACHE_RDS, "\n")
  } else {
    cat("[WARN] Cache file not written (no write permission or fetch failed)\n")
  }
} else {
  cat("[WARN] ClinVar data is empty. The ClinVar track will be blank in the app.\n")
  cat("       Check the messages above for the specific failure point.\n")
}

cat("\n")

# ── Step 3: Transcript exon structure cache ───────────────────────────────────

cat("--- Step 3: Transcript exon structure cache ---\n")

transcript_exons <- NULL

if (transcript_cache_is_fresh()) {
  message("Loading cached transcript exon structure from disk (< ",
          TRANSCRIPT_CACHE_DAYS, " days old)...")
  transcript_exons <- tryCatch(
    readRDS(TRANSCRIPT_CACHE_FILE),
    error = function(e) {
      cat("[FAIL] Could not read transcript cache file:", TRANSCRIPT_CACHE_FILE, "\n")
      cat("       Reason:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (!is.null(transcript_exons)) {
    cat("[OK]   Transcript cache loaded from disk:", nrow(transcript_exons), "exons\n")
  }
} else {
  message("Transcript cache absent or expired — fetching from Ensembl...")
  transcript_exons <- fetch_transcript_exons(TRANSCRIPT_ID)

  if (!is.null(transcript_exons)) {
    tryCatch({
      write_transcript_cache(transcript_exons)
      cat("[OK]   Transcript exon structure fetched and cached:",
          nrow(transcript_exons), "exons\n")
      cat("[OK]   Cache file written:", TRANSCRIPT_CACHE_FILE, "\n")
    }, error = function(e) {
      cat("[WARN] Exon structure fetched successfully but cache write failed.\n")
      cat("       Reason:", conditionMessage(e), "\n")
      cat("       This is likely a write-permission issue.\n")
    })
  } else {
    # Fetch failed — try stale cache as fallback
    if (file.exists(TRANSCRIPT_CACHE_FILE)) {
      warning(
        "Ensembl fetch failed. Falling back to stale transcript cache. ",
        "Coordinates may reflect an outdated annotation."
      )
      transcript_exons <- tryCatch(
        readRDS(TRANSCRIPT_CACHE_FILE),
        error = function(e) NULL
      )
      if (!is.null(transcript_exons)) {
        cat("[WARN] Ensembl fetch failed. Loaded stale cache as fallback.\n")
      } else {
        cat("[FAIL] Ensembl fetch failed and stale cache is unreadable.\n")
      }
    } else {
      cat("[FAIL] Ensembl fetch failed and no cache exists.\n")
      cat("       Coordinate mapping will produce NAs for all variants.\n")
    }
  }
}

cat("\n")

# ── Step 4: Summary ───────────────────────────────────────────────────────────

cat("================================================================================\n")
cat("  Summary\n")
cat("================================================================================\n")

clinvar_ok    <- !is.null(clinvar_data)    && nrow(clinvar_data) > 0
transcript_ok <- !is.null(transcript_exons) && nrow(transcript_exons) > 0

cat(sprintf("  ClinVar cache:    %s  (%s variants)\n",
    if (clinvar_ok)    "[OK]  " else "[FAIL]",
    if (clinvar_ok)    nrow(clinvar_data) else "0"))

cat(sprintf("  Transcript cache: %s  (%s exons)\n",
    if (transcript_ok) "[OK]  " else "[FAIL]",
    if (transcript_ok) nrow(transcript_exons) else "0"))

cat("\n")

if (clinvar_ok && transcript_ok) {
  cat("  Both caches are ready. The app should load without external API calls.\n")
} else if (!clinvar_ok && !transcript_ok) {
  cat("  Both fetches failed. Most likely cause: outbound HTTP is blocked on\n")
  cat("  this server. Confirm that eutils.ncbi.nlm.nih.gov and\n")
  cat("  rest.ensembl.org are reachable on port 443.\n")
} else if (!clinvar_ok) {
  cat("  ClinVar fetch failed. NCBI E-utilities may be blocked or rate-limited.\n")
  cat("  The ClinVar track will be empty. Internal variants should still display.\n")
} else {
  cat("  Transcript fetch failed. Ensembl REST API may be blocked.\n")
  cat("  Internal variants will display without genomic positions (all NA).\n")
}

cat("================================================================================\n\n")

################################################################################
# END OF CACHE BUILDER
################################################################################

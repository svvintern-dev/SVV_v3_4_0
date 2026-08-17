################################################################################
#                                                                              #
#                   SYNGAP1 VARIANT VIEWER - Shiny Application                 #
#                              Version 3.4.0                                   #
#                                                                              #
#  Purpose: Interactive genome browser for visualizing SYNGAP1 genetic         #
#           variants with patient data and research asset tracking             #
#                                                                              #
#  Architecture:                                                               #
#    1. Package Dependencies (~line 63)
#       - Package Installations
#       - Library Loading
#                                                                              #
#    2. ClinVar Data Loading (~line 125)                                       #
#       - Load ClinVar data from pre-built cache (see svv_build_cache.R)       #
#                                                                              #
#    3. Data Structures, Helpers & Coordinate Mapping (~lines 151-498)         #
#       - Amino acid code lookup table                                         #
#       - Protein nomenclature formatting functions                            #
#       - cDNA coordinate extraction                                           #
#       - Local cDNA→genome arithmetic (no per-coordinate API calls)           #
#       - Transcript exon structure loaded from pre-built cache                #
#         (see svv_build_cache.R)                                              #
#                                                                              #
#    4. Data Loading & Preprocessing (~line 500)                               #
#       - Load patient variant data from Endicott 2026_citizen_variants.csv    #       #
#       - Load 2026 Citizen registry from 2026_citizen_variants.xlsx           #
#         and aggregate per-variant patient counts                             #
#       - Load CureSyngap1 survey data from curesyngap1_variants.csv           #
#         and aggregate per-variant patient counts                             #
#       - Convert cDNA coordinates to genome positions via local arithmetic    #
#       - Merge positions and source counts into variant data frame            #
#                                                                              #
#    5. UI Definition (~line 779)                                              #
#       - Sidebar with dynamic variant type buttons                            #
#       - Research asset filter checkboxes (above ClinVar separator)           #
#       - Static ClinVar track button + log-spaced size filter slider (1–50,000 KB) #
#       - Main panel with IGV genome browser                                   #
#                                                                              #
#    6. GFF3 Track Formatting (~lines 974-1297)                                #
#       - Format internal variant data as GFF3 for IGV                         #
#       - Format ClinVar data as GFF3 for IGV                                  #
#                                                                              #
#    7. Server Logic (~line 1300)                                              #
#       - Reactive filtering based on user selections                          #
#       - Dynamic track loading/updating for internal variants                 #
#       - ClinVar track loading on button click                                #
#       - ClinVar track reloading on size slider change                        #
#       - IGV browser rendering and interaction                                #
#    
#    8. Run the app (~line 1845)                                                                        #
#                                                                              #
#  Design Philosophy:                                                          #
#    - Separate data processing from UI to improve performance                 #
#    - Cache external API calls to work offline and avoid rate limits          #
#    - Use reactive programming for smooth filter interactions                 #
#    - Generate tracks dynamically based on available variant types            #
#    - Expose genomic span (Size) in all feature popups                        #
#                                                                              #
################################################################################

################################################################################
# SECTION 1: PACKAGE DEPENDENCIES
################################################################################
#
#   Run svv_packages.R
#
#   That script installs all necessary packages for the app.
#   (This keeps app.R clean and avoids re-installing packages on every run.)
#
#===============================================================================
# PACKAGE INSTALLATION
#===============================================================================
#
#   This app depends on:
#   - Shiny (web framework)
#   - igvShiny (genome browser)
#   - dplyr/stringr (data manipulation)
#   - shinyWidgets (UI enhancements)
#   - readxl (Excel input)
#
# install.packages("shiny")
# install.packages("igvShiny")
# install.packages("dplyr")
# install.packages("stringr")
# install.packages("shinyWidgets")
# install.packages("readxl")
#
# If using Bioconductor packages:
# BiocManager::install("GenomicRanges")

#===============================================================================
# LIBRARY LOADING
#===============================================================================

# Core Shiny framework for web application
library(shiny)

# igvShiny: Wrapper for Integrative Genomics Viewer (IGV) in R
# Why: Provides interactive genome browser functionality directly in Shiny
library(igvShiny)

# Data manipulation with tidyverse syntax
# Why: Clean, readable syntax for filtering and transforming variant data
library(dplyr)

# String pattern matching and manipulation
# Why: Parse variant nomenclature (e.g., "c.333del", "p.R135X")
library(stringr)

# Extended Shiny widgets
# Why: sliderTextInput() provides a logarithmic-style slider by accepting an
#      explicit vector of stop values, giving fine-grained control at the low
#      end of the ClinVar size filter without sacrificing upper-range coverage
library(shinyWidgets)

# Read Excel files (.xlsx) for 2026 Citizen variant data
# Why: The 2026 Citizen registry is distributed as an Excel workbook;
#      readxl handles this natively without requiring Java (unlike xlsx package)
library(readxl)

################################################################################
# CITATION
################################################################################
# Shannon P, Gladki A, Scigocka K (2024). igvShiny: igvShiny: a wrapper of 
# Integrative Genomics Viewer (IGV - an interactive tool for visualization 
# and exploration integrated genomic data). R package version 1.2.0, 
# https://gladkia.github.io/igvShiny/, https://github.com/gladkia/igvShiny.
# ---------------------------------------------------------------------------
# END SECTION 1
# ---------------------------------------------------------------------------


################################################################################
# SECTION 2: CLINVAR DATA LOADING — FROM PRE-BUILT CACHE
################################################################################

# ClinVar data is loaded from a pre-built disk cache.
# Use svv_build_cache.R (run separately before launching the app) to fetch
# from NCBI and populate the cache. If the cache is absent, the ClinVar
# track will be empty this session.
#
# Cache files (written by svv_build_cache.R):
#   cache/clinvar/clinvar_data.rds   — parsed data frame
#   cache/clinvar/last_updated.txt   — ISO timestamp of last successful fetch

CLINVAR_CACHE_RDS <- "cache/clinvar/clinvar_data.rds"

clinvar_data <- if (file.exists(CLINVAR_CACHE_RDS)) {
  message("ClinVar: loading from cache...")
  readRDS(CLINVAR_CACHE_RDS)
} else {
  warning(
    "ClinVar cache not found at '", CLINVAR_CACHE_RDS, "'. ",
    "Run svv_build_cache.R to populate the cache. ",
    "The ClinVar track will be empty this session."
  )
  data.frame()
}
# ---------------------------------------------------------------------------
# END SECTION 2
# ---------------------------------------------------------------------------

################################################################################
# SECTION 3A: DATA STRUCTURES & CONSTANTS
################################################################################

# Amino Acid Code Mapping Table
# Purpose: Convert single-letter amino acid codes to three-letter codes
# Why: Three-letter codes are more readable for biologists
# Example: "R135X" → "Arg135ter" (Arginine at position 135 to stop codon)
aa_code_map <- c(
  A = "Ala",  # Alanine
  C = "Cys",  # Cysteine
  D = "Asp",  # Aspartate
  E = "Glu",  # Glutamate
  F = "Phe",  # Phenylalanine
  G = "Gly",  # Glycine
  H = "His",  # Histidine
  I = "Ile",  # Isoleucine
  K = "Lys",  # Lysine
  L = "Leu",  # Leucine
  M = "Met",  # Methionine
  N = "Asn",  # Asparagine
  P = "Pro",  # Proline
  Q = "Gln",  # Glutamine
  R = "Arg",  # Arginine
  S = "Ser",  # Serine
  T = "Thr",  # Threonine
  V = "Val",  # Valine
  W = "Trp",  # Tryptophan
  Y = "Tyr",  # Tyrosine
  X = "ter"   # Termination (stop codon)
)
# ---------------------------------------------------------------------------
# END SECTION 3A
# ---------------------------------------------------------------------------

################################################################################
# SECTION 3B: HELPER FUNCTIONS - PROTEIN NOMENCLATURE
################################################################################

#' Convert Single-Letter to Three-Letter Amino Acid Codes
#'
#' @description
#' Converts protein change notation from compact single-letter format to 
#' more readable three-letter format following HGVS nomenclature standards.
#'
#' @param protein_change_str String containing protein change (e.g., "p.R135X")
#'
#' @return String with three-letter amino acid codes (e.g., "Arg135ter")
#'
#' @details
#' Handles multiple protein change formats:
#'   - Simple substitutions: R135X → Arg135ter
#'   - Frameshifts: K138Hfs*11 → Lys138Hisfs*11
#'   - Synonymous: T120= → Thr120=
#'   - Deletions/insertions: maintained as-is
#'
#' Design Decision: Why do this conversion?
#'   - More readable for clinical researchers
#'   - Standard in clinical genetics publications
#'   - Reduces ambiguity (e.g., "Ser" vs "Cys" when using S vs C)
#'
#' @examples
#' swap_one_letter_to_three_letter("p.R135X")     # → "Arg135ter"
#' swap_one_letter_to_three_letter("R135X")       # → "Arg135ter" (p. prefix optional)
#' swap_one_letter_to_three_letter("K138Hfs*11")  # → "Lys138Hisfs*11"
swap_one_letter_to_three_letter <- function(protein_change_str) {
  # Remove the "p." prefix if present (HGVS standard prefix for protein changes)
  # Why: We want to work with just "R135X" not "p.R135X" for easier parsing
  protein_change_str <- gsub("^p\\.", "", protein_change_str)
  
  # Define regex pattern to match protein change components
  # Pattern breakdown:
  #   ^([A-Z])           - Start amino acid (single letter)
  #   ([0-9]+)           - Position number (e.g., 135)
  #   ([A-Z](?:fs\\*\\d+)|[A-Z]|=|del|ins|dup|[a-z]+\\*\\d+)$ - End change
  #
  # Why this pattern? Covers main variant types:
  #   - Simple changes: R135X (Arg→stop)
  #   - Frameshifts: K138Hfs*11
  #   - Synonymous: T120=
  #   - Special cases: del, ins, dup
  match <- str_match(protein_change_str, "^([A-Z])([0-9]+)([A-Z](?:fs\\*\\d+)|[A-Z]|=|del|ins|dup|[a-z]+\\*\\d+)$")
  
  # If pattern matches, proceed with conversion
  if (!is.na(match[1,1])) {
    # Extract components
    start_aa <- match[1,2]  # First amino acid (e.g., "R")
    position <- match[1,3]  # Position (e.g., "135")
    change <- match[1,4]    # Change (e.g., "X" for stop)
    
    # Convert start amino acid to three-letter code
    # Design: Use lookup table for fast, reliable conversion
    if (start_aa %in% names(aa_code_map)) {
      start_aa_three <- aa_code_map[start_aa]
    } else {
      # If not in map (shouldn't happen), keep original
      start_aa_three <- start_aa
    }
    
    # Handle the change component
    # Check if change is amino acid + frameshift/stop notation
    # Pattern: ([A-Z]|\*) for amino acid or stop, (fs\*\d+|\*\d+)? for frameshift
    change_match <- str_match(change, "^([A-Z]|\\*)(fs\\*\\d+|\\*\\d+)?$")
    
    if (!is.na(change_match[1,1])) {
      # Extract the amino acid and frameshift components
      change_aa <- change_match[1,2]    # Amino acid (e.g., "X" or "H")
      fs_part <- change_match[1,3]      # Frameshift notation (e.g., "fs*11")
      
      # Convert change amino acid to three-letter code
      if (change_aa %in% names(aa_code_map)) {
        change_aa_three <- aa_code_map[change_aa]
      } else {
        change_aa_three <- change_aa
      }
      
      # Reconstruct with frameshift notation if present
      # Example: "H" + "fs*11" → "Hisfs*11"
      change_three <- paste0(change_aa_three, ifelse(is.na(fs_part), "", fs_part))
    } else {
      # If change is 'fs*number' without preceding amino acid,
      # or other special notation like "=", "del", "ins"
      # Keep as-is since no amino acid to convert
      change_three <- change
    }
    
    # Reconstruct the full protein change string
    # Example: "Arg" + "135" + "ter" = "Arg135ter"
    protein_change_three <- paste0(start_aa_three, position, change_three)
    return(protein_change_three)
  } else {
    # If pattern doesn't match (e.g., complex variants, structural changes),
    # return original string unchanged
    # Design Decision: Fail gracefully rather than error out
    return(protein_change_str)
  }
}

#' Extract and Format Protein Change
#'
#' @description
#' Wrapper function that handles NA values and calls the conversion function.
#'
#' @param protein_change Raw protein change string from CSV
#'
#' @return Formatted protein change with three-letter codes, or NA if input is NA
#'
#' @details
#' This is a simple wrapper to handle edge cases before conversion.
#' Separates data validation from conversion logic (single responsibility principle).
extract_protein_change <- function(protein_change) {
  # Handle missing/empty values gracefully
  # Why: Not all variants have protein predictions (e.g., intronic variants)
  if (is.na(protein_change) || protein_change == "") return(protein_change)
  
  # Swap one-letter codes for three-letter codes
  protein_change_three <- swap_one_letter_to_three_letter(protein_change)
  return(protein_change_three)
}
# ---------------------------------------------------------------------------
# END SECTION 3B
# ---------------------------------------------------------------------------

################################################################################
# SECTION 3C: HELPER FUNCTIONS - COORDINATE EXTRACTION
################################################################################

#' Extract cDNA Coordinate from Variant String
#'
#' @description
#' Parses HGVS cDNA notation to extract the primary coordinate number.
#' Handles both single positions and ranges.
#'
#' @param variant HGVS cDNA notation (e.g., "c.333del", "c.190_200del")
#'
#' @return Numeric cDNA coordinate, or NA if unparseable
#'
#' @details
#' Handles multiple formats:
#'   - Simple: c.333del → 333
#'   - Range: c.190_200del → 190 (uses start position)
#'   - Intronic: c.190-2A>G → 190 (extracts exonic reference)
#'   - Invalid: c.GAIN → NA (no numeric coordinate)
#'
#' Design Decision: Why extract coordinates?
#'   - Needed for local cDNA→genome arithmetic (build_cdna_position_map)
#'   - cDNA coordinates are relative to transcript
#'   - IGV browser needs absolute genome coordinates
#'   - This bridges the gap between clinical nomenclature and visualization
#'
#' @examples
#' extract_cDNA_coordinate("c.333del")      # → 333
#' extract_cDNA_coordinate("c.190_200del")  # → 190 (start of range)
#' extract_cDNA_coordinate("c.190-2A>G")    # → 190 (intronic, uses exon ref)
#' extract_cDNA_coordinate("c.GAIN")        # → NA (no coordinate)
extract_cDNA_coordinate <- function(variant) {
  # Try to match simple format: c.NUMBER (e.g., c.333, c.490)
  # Pattern: c\. matches literal "c.", ([0-9]+) captures the number
  match <- str_match(variant, "c\\.([0-9]+)")
  
  if (!is.na(match[1, 2])) {
    # Successfully matched - return the coordinate as numeric
    # Why numeric? Needed for local coordinate mapping and sorting
    return(as.numeric(match[1, 2]))
  } else {
    # Didn't match simple format - try range format
    # Handle ranges like c.190_200 by taking the first part of the range
    # Why first position? 
    #   - Represents start of affected region
    #   - Local arithmetic maps single positions; we simplify ranges to a single point
    #   - IGV will show variant at this position (good enough for visualization)
    range_match <- str_match(variant, "c\\.([0-9]+)_([0-9]+)")
    
    if (!is.na(range_match[1, 2])) {
      # Return the first coordinate in the range as numeric
      return(as.numeric(range_match[1, 2]))
    } else {
      # Couldn't parse coordinate - return NA
      # This happens for:
      #   - Structural variants (e.g., "GAIN", "Entire coding sequence")
      #   - Complex rearrangements
      #   - Malformed strings
      # Design: Return NA rather than error - let downstream handle gracefully
      return(NA)
    }
  }
}
# ---------------------------------------------------------------------------
# END SECTION 3C
# ---------------------------------------------------------------------------

################################################################################
# SECTION 3D: LOCAL CDNA → GENOME COORDINATE MAPPING
################################################################################
#' Build a Complete cDNA-to-Genome Position Map Using Exon Structure
#'
#' @description
#' Given the exon data frame loaded from the transcript cache (see
#' svv_build_cache.R), walks across exons in transcript order and computes
#' in transcript order and computes the genomic coordinate for every cDNA
#' position in the dataset — entirely in R, with no further API calls.
#'
#' @param coords     Integer vector of cDNA coordinates to map (NAs are handled)
#' @param exon_df    Data frame of exon structure loaded from
#'                   cache/transcript_structure/ENST00000418600.rds
#'
#' @return Data frame with columns: cDNA_coordinate, seqid, start, end
#'         Unmappable coordinates (out of range or NA) get NA for start/end.
#'
#' @details
#' The mapping arithmetic:
#'
#'   cDNA position p falls in exon k when the cumulative transcript length
#'   through exon k-1 < p <= cumulative length through exon k.
#'
#'   Let offset = p - (cumulative length through exon k-1)
#'              = the 1-based position WITHIN exon k.
#'
#'   Plus-strand  (+1): genome position = exon_start + offset - 1
#'   Minus-strand (-1): genome position = exon_end   - offset + 1
#'
#'   Because Ensembl's /map/cdna/ endpoint returns a *range* of [start, end]
#'   for point queries (both equal the same base), we set start = end = the
#'   computed position, matching the previous behaviour exactly.
#'
#' Why c. coordinates start at 1, not 0:
#'   HGVS cDNA numbering is 1-based.  The coordinate extracted from "c.333del"
#'   is 333, meaning the 333rd coding base of the transcript.  This function
#'   treats all input coordinates as 1-based.
#'
#' Note on UTR / intronic variants:
#'   extract_cDNA_coordinate() already strips UTR prefixes (*) and intronic
#'   offsets (-N / +N), returning only the adjacent exonic coordinate.  Those
#'   simplified coordinates map correctly through this function.
#'
#' @examples
#' map   <- build_cdna_position_map(c(333, 490, NA, 99999), transcript_exons)
#' map   <- build_cdna_position_map(c(333, 490, NA, 99999), exons)
build_cdna_position_map <- function(coords, exon_df) {
  
  # Pre-compute cumulative transcript lengths so we can locate any coordinate
  # with a simple comparison rather than a loop-within-loop.
  #
  # cum_end[k]   = last cDNA position covered by exon k
  # cum_start[k] = first cDNA position covered by exon k
  cum_end   <- cumsum(exon_df$exon_length)
  cum_start <- c(1L, cum_end[-nrow(exon_df)] + 1L)
  strand    <- exon_df$strand[1]   # same for every row; grab once
  
  map_one <- function(p) {
    # Return a 2-element named integer for start and end genome position.
    # We keep start == end (point feature) to match the previous API behaviour.
    na_result <- c(start = NA_integer_, end = NA_integer_)
    
    if (is.na(p) || p < 1L) return(na_result)
    
    # Find which exon contains cDNA position p
    hit <- which(p >= cum_start & p <= cum_end)
    
    if (length(hit) == 0L) {
      # Coordinate exceeds transcript length — out of range
      message("cDNA coordinate ", p, " is out of transcript range (max ",
              max(cum_end), "); skipping.")
      return(na_result)
    }
    
    k      <- hit[1]                          # exon index (1-based)
    offset <- p - cum_start[k] + 1L          # 1-based position within this exon
    
    if (strand == 1L) {
      # Plus strand: cDNA increases in the same direction as genomic coordinates
      gpos <- exon_df$exon_start[k] + offset - 1L
    } else {
      # Minus strand: cDNA increases as genomic coordinates decrease
      gpos <- exon_df$exon_end[k] - offset + 1L
    }
    
    return(c(start = gpos, end = gpos))
  }
  
  # Apply to all unique coordinates (vectorised over the input vector)
  unique_coords <- unique(coords)
  mapped        <- lapply(unique_coords, map_one)
  
  data.frame(
    cDNA_coordinate = unique_coords,
    seqid           = "chr6",
    start           = sapply(mapped, `[[`, "start"),
    end             = sapply(mapped, `[[`, "end"),
    stringsAsFactors = FALSE
  )
}


# ---------------------------------------------------------------------------
# Load transcript exon structure from pre-built cache.
# Use svv_build_cache.R (run separately before launching the app) to fetch
# from Ensembl and populate the cache. If the cache is absent, all
# cDNA→genome coordinate mappings will be NA this session.
#
# Cache files (written by svv_build_cache.R):
#   cache/transcript_structure/ENST00000418600.rds           — exon data frame
#   cache/transcript_structure/ENST00000418600.rds.timestamp — ISO fetch time
# ---------------------------------------------------------------------------

TRANSCRIPT_CACHE_FILE <- "cache/transcript_structure/ENST00000418600.rds"

transcript_exons <- if (file.exists(TRANSCRIPT_CACHE_FILE)) {
  message("Loading transcript exon structure from cache...")
  readRDS(TRANSCRIPT_CACHE_FILE)
} else {
  warning(
    "Transcript cache not found at '", TRANSCRIPT_CACHE_FILE, "'. ",
    "Run svv_build_cache.R to populate the cache. ",
    "All cDNA→genome coordinate mappings will be NA this session."
  )
  NULL
}
# ---------------------------------------------------------------------------
# END SECTION 3D
# ---------------------------------------------------------------------------

################################################################################
# SECTION 4: DATA LOADING & PREPROCESSING
################################################################################

# ------------------------------------------------------------------------------
#  Load Variant Data
# ------------------------------------------------------------------------------
#' @description
#' Loads patient variant data from CSV file.
#'
#' Design Decision: Why load at startup rather than in server?
#'   - Data doesn't change during app session
#'   - Loading once is more efficient than per-user
#'   - Makes debugging easier (data loaded before UI interactions)
#'   - Preprocessing happens once, not per user connection
#'
#' File Format: Endicott 2026_citizen_variants.csv
#'   - 153 rows (variants)
#'   - 33 columns (metadata, patient info, research assets)
#'   - Key columns:
#'     - SYNGAP1.variant: cDNA notation (e.g., "c.333del")
#'     - variant.type: Variant classification
#'     - predicted.protein: Protein change
#'     - biorepository, iPSC.line, mouse.line: Research assets
#'     - X..patients.in.Citizen.Health: 2024 Citizen patient count

# Load CSV data
# stringsAsFactors = FALSE prevents automatic factor conversion
# Why? We want to manipulate strings freely without factor level constraints
variant_data <- read.csv("Endicott 2026_citizen_variants.csv", stringsAsFactors = FALSE)


curesyngap1_counts <- tryCatch({
  raw <- read.csv("curesyngap1_variants.csv", stringsAsFactors = FALSE)
  table(raw$variant)
}, error = function(e) {
  warning(
    "CureSyngap1 data could not be loaded from 'curesyngap1_variants.csv': ",
    conditionMessage(e), ". ",
    "CureSyngap1 Survey Count will show as 0 for all variants this session."
  )
  table(character(0))
})
# ==============================================================================
# Standardize Variant Type Names
# ==============================================================================
#' @description
#' Normalizes variant type nomenclature for consistency.
#'
#' Problem: CSV has inconsistent naming
#'   - "missense VUS" vs "VUS - missense"
#'   - "frameshift" vs "frameshift deletion" vs "frameshift insertion"
#'   - Mixed case
#'
#' Solution: Map all variations to canonical names
#'
#' Why important?
#'   - Consistent button labels in UI
#'   - Reliable filtering
#'   - Correct color assignments
#'   - Clean track names in IGV

variant_data <- variant_data %>%
  mutate(
    # First, convert all to lowercase for case-insensitive matching
    variant.type = tolower(variant.type),
    
    # Then map to canonical names using case_when (like switch statement)
    # Order matters: more specific patterns first
    variant.type = case_when(
      # Missense variants (amino acid substitutions)
      variant.type %in% c('missense') ~ 'missense',
      
      # Nonsense variants (premature stop codons)
      variant.type %in% c('nonsense') ~ 'nonsense',
      
      # Frameshift variants (reading frame disruptions)
      # Groups all frameshift subtypes together
      variant.type %in% c('frameshift', 'frameshift deletion', 'frameshift insertion') ~ 'frameshift',
      
      # Splice (donor or acceptor)
      variant.type %in% c('splice','splice acceptor','splice donor') ~ 'splice',
      
      # Indels (small insertions/deletions maintaining frame)
      variant.type %in% c('indel', 'insertion deletion', 'in-frame','in-frame delins','in-frame deletion') ~ 'indel',
      
      # CNV (large multigenic deletions and single exons removed)
      # CNV variants are currently hidden because spread sheet does not contain CNV locations
      # this can be updated by removing comments from the cnv button in two parts of section 7
      variant.type %in% c('cnv') ~ 'cnv',
      
      # VUS (general uncertain significance)
      variant.type %in% c('vus','missense vus','missense-vus', 'nonsense vus') ~ 'vus',
      
      # Other variants (mostly intronic)
      variant.type %in% c('intronic','in-frame duplication','synonymous') ~ 'other',
      
      # Default: catch all for things without match
      # Handles new types without code changes
      TRUE ~ 'other'
    )
  )

# Convert to factor for efficient storage and categorical operations
# Why factor? 
#   - Memory efficient (stores as integers with labels)
#   - Enforces valid values (can't accidentally add typos)
#   - Useful for grouping and counting
variant_data$variant.type <- factor(variant_data$variant.type)

# ==============================================================================
#  Precompute Derived Fields
# ==============================================================================
#' @description
#' Extract and format data needed for visualization before app starts.
#'
#' Design Philosophy: Do expensive work once at startup
#'   - Parsing happens once, not every time a user filters
#'   - Results stored in data frame for fast access
#'   - Users get instant UI response
#'
#' What's being precomputed:
#'   1. cDNA coordinates (for local exon-structure mapping)
#'   2. Formatted protein changes (for display)

variant_data <- variant_data %>%
  mutate(
    # Extract numeric cDNA coordinate from variant string
    # sapply applies function to each row (vectorized operation)
    # Example: "c.333del" → 333
    cDNA_coordinate = sapply(SYNGAP1.variant, extract_cDNA_coordinate),
    
    # Format protein changes with three-letter amino acids
    # Example: "p.R135X" → "Arg135ter"
    protein_change = sapply(predicted.protein, extract_protein_change),
    
    # Per-variant patient count from the CureSyngap1 survey export.
    # Same lookup pattern as above.
    curesyngap1_count = as.integer(
      ifelse(
        SYNGAP1.variant %in% names(curesyngap1_counts),
        curesyngap1_counts[SYNGAP1.variant],
        0L
      )
    )
  )

#' Map cDNA Coordinates to Genome Positions (Local Arithmetic)
#'
#' @description
#' Converts all cDNA coordinates to genome positions using the pre-fetched
#' transcript exon structure.  No further API calls are made here — the
#' mapping is pure arithmetic over the exon table.
#'
#' Performance:
#'   - Reads cached exon .rds from disk (~1 ms) + local arithmetic (~1 ms)
#'   - If cache is absent: run svv_build_cache.R to fetch from Ensembl
#'   - Compare to old approach: 143 API calls × ~300 ms = ~43 s on cold cache

if (!is.null(transcript_exons)) {
  # Local mapping: no API calls, typically completes in < 1 ms
  positions_df <- build_cdna_position_map(
    coords   = variant_data$cDNA_coordinate,
    exon_df  = transcript_exons
  )
} else {
  # Fallback: fill with NAs so the rest of the pipeline runs unchanged
  warning(
    "Transcript exon structure unavailable. ",
    "Genome positions will be NA — variants will not display in IGV."
  )
  unique_coords <- unique(variant_data$cDNA_coordinate)
  positions_df  <- data.frame(
    cDNA_coordinate = unique_coords,
    seqid           = "chr6",
    start           = NA_integer_,
    end             = NA_integer_,
    stringsAsFactors = FALSE
  )
}

# Merge genome positions back into main variant data.
# Left join keeps all variants; those without a valid position get NA.
# Result: variant_data gains seqid, start, end columns for IGV display.
variant_data <- merge(variant_data, positions_df, by = "cDNA_coordinate")

#' Define Color Scheme for Variant Types
#'
#' @description
#' Color palette for different mutation types in IGV browser.
#'
#' Design Considerations:
#'   - Colorblind-friendly (avoid red/green for critical distinctions)
#'   - Sufficient contrast for visibility
#'   - Semantic meaning where possible
#'   - Consistent with scientific conventions
#'
#' Color Choices:
#'   - Blue (#3E6CB5): Missense - common, moderate impact
#'   - Orange (#FFA500): Missense-VUS - uncertain, needs attention
#'   - Green (#19AE69): Nonsense - clearly pathogenic
#'   - Purple (#6F3592): Frameshift - severe disruption
#'   - Light Blue (#D0DDEE): Indel - variable impact
#'   - Light Green (#CBE7D4): Gain - structural variant
#'   - Lavender (#D4CAE1): VUS - uncertain
#'   - Gray (#4E4E4E): Intronic/Loss - likely benign or uncertain

color_table <- list(
  missense   = "#3E6CB5",      # Blue
  nonsense   = "#19AE69",      # Green
  frameshift = "#6F3592",      # Purple
  splice     = "#D4CAE1",      # Lavender
  indel      = "#D0DDEE",      # Light Blue
  cnv        = "#CBE7D4",      # Light Green
  vus        = "#FFA500",      # Orange (attention-grabbing for VUS)
  other      = "#4E4E4E",      # Gray
  clinvar    = "#00BCD4"       # Teal - visually distinct from all internal track colors
)

# ---------------------------------------------------------------------------
# END SECTION 4
# ---------------------------------------------------------------------------

################################################################################
# SECTION 5: USER INTERFACE DEFINITION
################################################################################

#' Shiny UI Layout
#'
#' @description
#' Defines the visual layout and interactive elements of the web application.
#'
#' Architecture: Two-column layout
#'   - Left sidebar (30% width): Controls and filters
#'   - Right main panel (70% width): IGV genome browser
#'
#' Design Philosophy:
#'   - Put controls on left (common web convention)
#'   - Maximize space for visualization
#'   - Group related controls together
#'   - Clear visual hierarchy

ui <- shinyUI(fluidPage(
  # Application title
  # Appears at top of page
  titlePanel("SVV for SYNGAP1 Variant Viewer"),
  
  # Two-column layout with sidebar and main content
  sidebarLayout(
    
    # ============================================================
    # LEFT SIDEBAR: Controls
    # ============================================================
    sidebarPanel(
      # Section header for variant type buttons
      h3("Add Tracks"),
      
      # Dynamic UI element: Buttons generated based on available variant types
      # Why dynamic? 
      #   - Number of variant types varies by dataset
      #   - Filters change which types are visible
      #   - Automatically adapts to data
      # Rendered by server (see output$trackButtons below)
      uiOutput("trackButtons"),
      
      # ========================================
      # Research Asset Filters
      # ========================================
      # Placed above the ClinVar separator so all internal-data controls
      # (variant type buttons + these checkboxes) are grouped together,
      # with the external ClinVar section clearly separated below.
      #
      # Why these filters?
      #   - Help researchers find variants with available resources
      #   - Enables focused analysis on well-characterized variants
      #   - Supports research planning and collaboration
      
      # Filter 1: Biorepository samples
      # Shows only variants where biological samples are stored
      # Use case: "Which variants can we get DNA/RNA for?"
      div(
        class = "form-group",
        checkboxInput("biorepository_filter", 
                      "Has biorepository samples", 
                      value = FALSE)  # Default unchecked (show all)
      ),
      
      # Filter 2: iPSC cell lines
      # Shows only variants with induced pluripotent stem cell lines
      # Use case: "Which variants can we study in cell culture?"
      div(
        class = "form-group",
        checkboxInput("iPSC_line_filter", 
                      "Cell line available", 
                      value = FALSE)
      ),
      
      # Filter 3: Mouse models
      # Shows only variants with corresponding mouse lines
      # Use case: "Which variants have animal models for testing therapies?"
      div(
        class = "form-group",
        checkboxInput("mouse_line_filter", 
                      "Mouse line available", 
                      value = FALSE)
      ),
      
      # ========================================
      # ClinVar Track Button (static)
      # ========================================
      # This button is static (not part of the dynamic variant-type loop)
      # because ClinVar is an external reference dataset, not a filtered
      # subset of the internal Citizen Health variant data.
      # It always appears regardless of which filters are active.
      hr(style = "border-top: 1px solid #ccc; margin: 8px 0;"),
      actionButton(
        inputId = "addTrack_clinvar",
        label   = "ClinVar"
      ),
      
      # ClinVar size filter slider
      # Filters ClinVar variants shown on the track by their genomic span.
      #
      # Why sliderTextInput() instead of sliderInput()?
      #   ClinVar spans an enormous size range: ~1 bp SNVs up to ~46,632 KB
      #   whole-chromosome CNVs. A linear slider across that range would make
      #   the low end (where most interesting variants live) nearly impossible
      #   to control. sliderTextInput() accepts an explicit vector of stop
      #   values, so we can space them logarithmically: fine steps at the bottom
      #   (1, 2, 3 … 10 KB) and progressively coarser steps toward the top
      #   (5,000 KB steps above 10,000 KB). Every stop is a "round" KB value
      #   that makes clinical sense.
      #
      # Stop sequence (45 positions total):
      #   1–9 KB      : step 1 KB   (fine control for point variants / small indels)
      #   10–90 KB    : step 10 KB  (single-exon to multi-exon deletions)
      #   100–900 KB  : step 100 KB (sub-megabase SVs)
      #   1,000–9,000 KB : step 1,000 KB (megabase-scale CNVs)
      #   10,000–50,000 KB : step 5,000 KB (cytogenetic-band CNVs)
      #
      # Default: 100 KB — excludes whole-gene CNVs while retaining variants
      #   within the SYNGAP1 locus (~37 kb). Users drag right to reveal
      #   progressively larger structural variants.
      #
      # Behavior: If the ClinVar track is already displayed, moving the slider
      # immediately reloads it with only variants at or below the new threshold.
      sliderTextInput(
        inputId  = "clinvar_size_kb",
        label    = "Max ClinVar variant size (KB)",
        choices  = c(
          seq(1,     9,     by = 1),      # 1–9 KB     : step 1 KB
          seq(10,    90,    by = 10),     # 10–90 KB   : step 10 KB
          seq(100,   900,   by = 100),    # 100–900 KB : step 100 KB
          seq(1000,  9000,  by = 1000),   # 1–9 MB     : step 1,000 KB
          seq(10000, 50000, by = 5000)    # 10–50 MB   : step 5,000 KB
        ),
        selected = 100,       # Default: 100 KB
        grid     = FALSE,
        width    = "100%"
      ),
      
      # Set sidebar width to 30% of page (3 out of 12 columns)
      # Why 30%? 
      #   - Enough space for controls without crowding
      #   - Leaves 70% for genome browser (where the action is)
      width = 3
    ),
    
    # ============================================================
    # RIGHT MAIN PANEL: Visualization
    # ============================================================
    mainPanel(
      # IGV genome browser widget
      # This is where variants are visualized on chromosome 6
      #
      # Parameters:
      #   - id: 'igvShiny_0' - unique identifier for this instance
      #   - height: "700px" - tall enough to see multiple tracks
      #
      # Rendered by server (see output$igvShiny_0 below)
      igvShinyOutput('igvShiny_0', height = "700px"),
      
      # Set main panel width to 70% of page (9 out of 12 columns)
      width = 9
    )
  ),
  
  # ============================================================
  # FOOTER: Citation
  # ============================================================
  # Horizontal rule for visual separation
  hr(),
  
  # Citation div with centered text and small font
  # Why include citation?
  #   - Give credit to igvShiny developers
  #   - Help users find original tool for their own projects
  #   - Professional appearance
  tags$div(
    style = "text-align: center; font-size: 12px;",
    tags$p("Citation:"),
    tags$p(
      "Shannon P, Gladki A, Scigocka K (2024). ",
      tags$em("igvShiny: igvShiny: a wrapper of Integrative Genomics Viewer "),
      "(IGV - an interactive tool for visualization and exploration integrated genomic data). ",
      "R package version 1.2.0, ",
      tags$a(href = "https://gladkia.github.io/igvShiny/", 
             "https://gladkia.github.io/igvShiny/"), ", ",
      tags$a(href = "https://github.com/gladkia/igvShiny", 
             "https://github.com/gladkia/igvShiny.")
    )
  )
))
# ---------------------------------------------------------------------------
# END SECTION 5
# ---------------------------------------------------------------------------

################################################################################
# SECTION 6A: DATA FORMATTING FOR IGV
################################################################################

#' Create GFF3-formatted Data for IGV Tracks
#'
#' @description
#' Converts variant data to GFF3 format required by IGV browser.
#'
#' @param data Filtered variant data frame
#' @param variant_type Type of variant for this track (e.g., "missense")
#'
#' @return Data frame in GFF3 format for igvShiny
#'
#' @details
#' GFF3 Format Requirements:
#'   1. seqid    - Chromosome (chr6)
#'   2. source   - Data source (not used, set to ".")
#'   3. type     - Feature type (variant type)
#'   4. start    - Genomic start position
#'   5. end      - Genomic end position
#'   6. score    - Confidence score (not used, set to ".")
#'   7. strand   - DNA strand (+/-)
#'   8. phase    - Coding phase (not used for variants, set to ".")
#'   9. attributes - Semicolon-separated key=value pairs
#'
#' Attributes Field:
#'   Contains metadata displayed when user clicks variant:
#'   - ID: Unique identifier
#'   - Name: Display name (protein change if available, else cDNA)
#'   - Description: Full cDNA notation
#'   - Number of patients: Clinical relevance indicator
#'   - Size: Genomic span (bp or kb)
#'
#' Design Decision: Why GFF3?
#'   - Standard genomics format (widely understood)
#'   - Supported natively by IGV
#'   - Flexible attributes system
#'   - Human-readable text format

create_gff3_data <- function(data, variant_type) {
  data <- data %>%
    # Remove variants without genome positions
    # Why filter? IGV can't display variants without coordinates
    # These are typically structural variants or failed API queries
    filter(!is.na(start)) %>%
    
    mutate(
      # Column 1: Chromosome
      seqid = "chr6",
      
      # Column 2: Source
      # "." indicates no specific source
      # Could use "patient_database" or "clinical_variants" if needed
      source = ".",
      
      # Column 3: Feature type
      # Uses the variant type (missense, nonsense, etc.)
      # This is what gets colored in IGV based on color_table
      type = variant_type,
      
      # Column 6: Score
      # "." indicates no score
      # Could add CADD scores, conservation scores, etc. if available
      score = ".",
      
      # Column 7: Strand
      # "+" is used here as a GFF3 display convention for IGV feature rendering.
      # Note: SYNGAP1 is actually transcribed on the minus (−) strand of chr6;
      # however, the genomic start/end coordinates we supply are already
      # absolute hg38 positions (computed correctly accounting for strand in
      # build_cdna_position_map), so the strand field here does not affect
      # coordinate accuracy — it only influences how IGV draws arrow direction.
      strand = "+",
      
      # Column 8: Phase
      # "." for not applicable (phase is for CDS features, not variants)
      phase = ".",
      
      # Column 9: Attributes (most important for display!)
      # Build different attribute strings based on data availability
      attributes = dplyr::if_else(
        # Check if protein change is available
        is.na(protein_change) | protein_change == "",
        
        # NO protein change - use cDNA notation only
        # This happens for intronic variants, synonymous variants, etc.
        paste0(
          "ID=", SYNGAP1.variant,           # Unique ID
          ";Name=", SYNGAP1.variant,        # Display name
          ";Description=cDNA: ", SYNGAP1.variant,  # Tooltip text
          ";2026CitizenCount=", X..patients.in.Citizen.Health,
          ";CURESYNGAP1SurveyCount=", curesyngap1_count,
          ";type=", variant_type,
          ";Size=", dplyr::if_else(                 # Genomic span in bp (or kb if >= 1000)
            pmax(end - start + 1, 1L) >= 1000,
            paste0(round(pmax(end - start + 1, 1L) / 1000, 2), " kb"),
            paste0(pmax(end - start + 1, 1L), " bp")
          )
        ),
        
        # YES protein change - show protein change as name
        # This is what users typically want to see
        paste0(
          "ID=", protein_change,            # Unique ID (protein level)
          ";Name=", protein_change,         # Display name (e.g., "Arg135ter")
          ";Description=cDNA: ", SYNGAP1.variant,  # Still show cDNA in tooltip
          ";2026CitizenCount=", X..patients.in.Citizen.Health,          
          ";CURESYNGAP1SurveyCount=", curesyngap1_count,
          ";type=", variant_type,
          ";Size=", dplyr::if_else(                 # Genomic span in bp (or kb if >= 1000)
            pmax(end - start + 1, 1L) >= 1000,
            paste0(round(pmax(end - start + 1, 1L) / 1000, 2), " kb"),
            paste0(pmax(end - start + 1, 1L), " bp")
          )
        )
      )
    ) %>%
    # Select and order columns to match GFF3 specification
    # Order is critical - IGV expects columns in this exact order
    select(seqid, source, type, start, end, score, strand, phase, attributes)
  
  return(data)
}
# ---------------------------------------------------------------------------
# END SECTION 6A
# ---------------------------------------------------------------------------

################################################################################
# SECTION 6B: CLINVAR DATA FORMATTING FOR IGV
################################################################################

#' Create GFF3-formatted Data from ClinVar API Data for IGV Tracks
#'
#' @description
#' Converts the parsed ClinVar data frame into GFF3 format required by IGV.
#' All available ClinVar metadata fields are packed into the GFF3 attributes
#' column so they appear in IGV's feature detail pop-up when a user clicks
#' a variant.
#'
#' @param data The clinvar_data data frame (fetched at startup via NCBI E-utilities)
#' @param max_size_kb Numeric. Maximum variant span in kilobases to include in
#'   the track. Variants whose (end - start + 1) exceeds this threshold are excluded.
#'   Defaults to 100 KB. Driven by the "Max ClinVar variant size" slider in the
#'   UI; pass as.numeric(input$clinvar_size_kb) from the server.
#'
#' @return Data frame in GFF3 format for use with loadGFF3TrackFromLocalData()
#'
#' @details
#' GFF3 columns produced:
#'   seqid      - "chr6" (all SYNGAP1 variants are on chromosome 6)
#'   source     - "ClinVar"
#'   type       - "clinvar" (drives teal color via color_table)
#'   start      - GRCh38 start position (parsed from GRCh38Location)
#'   end        - GRCh38 end position
#'   score      - "."
#'   strand     - "+"
#'   phase      - "."
#'   attributes - All ClinVar metadata fields as key=value pairs
#'
#' Attributes included (all available ClinVar fields):
#'   Name, ProteinChange, Condition, Accession, VariationID, AlleleID,
#'   dbSNP_ID, CanonicalSPDI, VariantType, MolecularConsequence,
#'   GermlineClassification, GermlineLastEvaluated, GermlineReviewStatus,
#'   SomaticClinicalImpact, SomaticLastEvaluated, SomaticReviewStatus,
#'   OncogenicityClassification, OncogenicityLastEvaluated,
#'   OncogenicityReviewStatus, Size
#'
#' Design Decision: Why include all fields?
#'   - ClinVar metadata is rich and clinically meaningful
#'   - Users can click a variant to see full clinical context
#'   - Avoids information loss from the source dataset
#'   - Mirrors IGV's standard ClinVar track behavior

create_clinvar_gff3_data <- function(data, max_size_kb = 100) {
  
  # Guard: return empty GFF3 frame if ClinVar data failed to load.
  # This happens when both the NCBI fetch and the cache are unavailable,
  # leaving clinvar_data as a zero-column data.frame().
  # Without this guard, dplyr::filter() crashes trying to reference
  # clinvar_start on a frame that has no columns at all.
  required_cols <- c("clinvar_start", "clinvar_end", "clinvar_chr")
  if (nrow(data) == 0 || !all(required_cols %in% names(data))) {
    message("ClinVar: no data available — track will be empty.")
    return(data.frame())
  }
  
  # Helper: sanitize a value for GFF3 attribute string
  # Replaces semicolons and equals signs which are GFF3 delimiters,
  # and trims whitespace. Returns "." for NA/empty values.
  clean_attr <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x[is.na(x) | x == ""] <- "."
    # Escape GFF3 reserved characters inside values
    x <- gsub(";", "%3B", x, fixed = TRUE)
    x <- gsub("=", "%3D", x, fixed = TRUE)
    x
  }
  
  # Convert KB threshold to base pairs for coordinate arithmetic
  max_size_bp <- max_size_kb * 1000
  
  data %>%
    # Remove variants without valid GRCh38 coordinates
    # These are typically older submissions that predate hg38 mapping
    filter(!is.na(clinvar_start) & !is.na(clinvar_end)) %>%
    
    # Remove variants whose span exceeds the user-specified size threshold.
    #
    # Why: Large copy number variants (e.g. chr6:156,974-46,789,291) span
    # millions of base pairs and physically cover the entire SYNGAP1 locus.
    # At every zoom level they sit on top of all smaller point variants,
    # intercepting every click and making the underlying variants unreachable.
    #
    # The threshold is controlled by the "Max variant size" slider in the UI
    # (default 100 kb). Lowering it progressively hides larger SVs and CNVs,
    # letting the user focus on point variants and small indels.
    #   - SYNGAP1 itself is ~37 kb, so 100 kb already excludes whole-gene CNVs.
    #   - Large CNVs spanning cytogenetic bands are not interpretable at this
    #     zoom level and are better reviewed directly in ClinVar's own browser.
    filter((clinvar_end - clinvar_start) <= max_size_bp) %>%
    
    mutate(
      # --- Required GFF3 columns ---
      seqid  = clinvar_chr,          # Chromosome (e.g. "chr6")
      source = "ClinVar",            # Data provenance label
      type   = "clinvar",            # Drives teal color via color_table
      start  = clinvar_start,
      end    = clinvar_end,
      score  = ".",
      strand = "+",                  # GFF3 display convention (see note in create_gff3_data)
      phase  = ".",
      
      # --- Attributes column: all ClinVar metadata ---
      # IGV displays these as a table when the user clicks a feature.
      # Fields are semicolon-separated key=value pairs per GFF3 spec.
      attributes = paste0(
        "ID=",                              clean_attr(Accession),
        ";Name=",                           clean_attr(Name),
        ";ProteinChange=",                  clean_attr(Protein.change),
        ";Condition=",                      clean_attr(Condition.s.),
        ";Accession=",                      clean_attr(Accession),
        ";VariationID=",                    clean_attr(VariationID),
        ";AlleleID=",                       clean_attr(AlleleID.s.),
        ";dbSNP_ID=",                       clean_attr(dbSNP.ID),
        ";CanonicalSPDI=",                  clean_attr(Canonical.SPDI),
        ";VariantType=",                    clean_attr(Variant.type),
        ";MolecularConsequence=",           clean_attr(Molecular.consequence),
        ";GermlineClassification=",         clean_attr(Germline.classification),
        ";GermlineLastEvaluated=",          clean_attr(Germline.date.last.evaluated),
        ";GermlineReviewStatus=",           clean_attr(Germline.review.status),
        ";SomaticClinicalImpact=",          clean_attr(Somatic.clinical.impact),
        ";SomaticLastEvaluated=",           clean_attr(Somatic.clinical.impact.date.last.evaluated),
        ";SomaticReviewStatus=",            clean_attr(Somatic.clinical.impact.review.status),
        ";OncogenicityClassification=",     clean_attr(Oncogenicity.classification),
        ";OncogenicityLastEvaluated=",      clean_attr(Oncogenicity.date.last.evaluated),
        ";OncogenicityReviewStatus=",       clean_attr(Oncogenicity.review.status),
        ";Size=",                           ifelse(
          pmax(clinvar_end - clinvar_start, 1L) >= 1e6,
          paste0(round(pmax(clinvar_end - clinvar_start, 1L) / 1e6, 2), " Mb"),
          ifelse(
            pmax(clinvar_end - clinvar_start, 1L) >= 1000,
            paste0(round(pmax(clinvar_end - clinvar_start, 1L) / 1000, 2), " kb"),
            paste0(pmax(clinvar_end - clinvar_start, 1L), " bp")
          )
        )
      )
    ) %>%
    # Select GFF3 columns in the required order
    select(seqid, source, type, start, end, score, strand, phase, attributes)
}
# ---------------------------------------------------------------------------
# END SECTION 6B
# ---------------------------------------------------------------------------

################################################################################
# SECTION 6C: CLINVAR TRACK PINNING HELPER
################################################################################

#' Re-pin the ClinVar Track to the Bottom of the IGV Stack
#'
#' @description
#' Removes the ClinVar track from IGV and immediately re-adds it so that it
#' always sits below every other variant track, regardless of add order.
#'
#' Why needed?
#'   IGV appends new tracks to the bottom of its track list.  If ClinVar is
#'   already displayed when the user clicks a variant-type button, the new
#'   track would appear below ClinVar — pushing it above the new track.
#'   Removing and re-loading ClinVar forces it to the bottom of the stack.
#'
#' @param session  Shiny session object
#' @param data     The clinvar_data data frame (fetched at startup)
#' @param size_kb  Current slider value (character or numeric KB threshold)

repin_clinvar_to_bottom <- function(session, data, size_kb) {
  # Remove existing ClinVar track (if present) so IGV forgets its position
  removeTracksByName(session, id = "igvShiny_0", "ClinVar")
  
  # Re-build GFF3 with the current size threshold
  clinvar_gff3 <- create_clinvar_gff3_data(data,
                                           max_size_kb = as.numeric(size_kb))
  
  if (nrow(clinvar_gff3) > 0) {
    loadGFF3TrackFromLocalData(
      session,
      id               = "igvShiny_0",
      trackName        = "ClinVar",
      tbl              = clinvar_gff3,
      colorTable       = color_table,
      colorByAttribute = "type",
      trackHeight      = 50,
      displayMode      = "EXPANDED",
      visibilityWindow = 1e8
    )
  } else {
    # All ClinVar variants exceed the threshold — re-add empty track
    loadGFF3TrackFromLocalData(
      session,
      id               = "igvShiny_0",
      trackName        = "ClinVar",
      tbl              = data.frame(),
      colorTable       = color_table,
      colorByAttribute = "type",
      trackHeight      = 50,
      displayMode      = "EXPANDED",
      visibilityWindow = 1e8
    )
  }
}
# ---------------------------------------------------------------------------
# END SECTION 6C
# ---------------------------------------------------------------------------

################################################################################
# SECTION 7: SERVER LOGIC (REACTIVE PROGRAMMING)
################################################################################

#' Shiny Server Function
#'
#' @description
#' Defines the reactive logic that powers the application.
#'
#' Reactive Programming Concepts:
#'   - Reactive values: Change over time (e.g., checkbox state)
#'   - Reactive expressions: Auto-recalculate when inputs change
#'   - Observers: Execute side effects when inputs change
#'
#' Server Architecture:
#'   1. State management (track which buttons were clicked)
#'   2. Data filtering (respond to checkbox changes)
#'   3. UI generation (create buttons based on filtered data)
#'   4. Event handling (respond to button clicks)
#'   5. Track updates (refresh IGV when filters change)
#'   6. IGV initialization (set up genome browser)

server <- function(input, output, session) {
  
  # ============================================================
  # STATE MANAGEMENT
  # ============================================================
  
  #' Track State Tracker
  #'
  # @description
  # Keeps track of which variant type tracks have been added to IGV.
  #
  # Why needed?
  #   - Remember user's track selections across filter changes
  #   - Avoid re-adding tracks that are already displayed
  #   - Enable smart updates (only refresh visible tracks)
  #
  # Data structure: Named list
  #   Example: list(missense = TRUE, nonsense = TRUE, frameshift = FALSE)
  #   TRUE = track is currently displayed
  #   FALSE or missing = track not displayed
  #
  # reactiveVal() creates a reactive variable that can be read and set
  # Like a reactive version of a regular variable
  added_tracks <- reactiveVal(list())
  
  # ClinVar Track State Tracker
  #
  # Separate from added_tracks because ClinVar is not part of the dynamic
  # variant-type loop. This flag lets the slider observer know whether to
  # reload the ClinVar track when the size threshold changes.
  #   TRUE  = ClinVar track is currently displayed in IGV
  #   FALSE = ClinVar track has not been added yet (slider changes are no-ops)
  clinvar_track_added <- reactiveVal(FALSE)
  
  # ============================================================
  # DATA FILTERING (REACTIVE EXPRESSION)
  # ============================================================
  
  #' Filtered Variant Data
  #'
  # @description
  # Reactive expression that filters variant data based on checkbox selections.
  #
  # Reactive Expression Behavior:
  #   - Automatically recalculates when inputs change
  #   - Caches result until inputs change (efficient)
  #   - Can be used by multiple observers/outputs
  #
  # Why reactive()? 
  #   - Avoid duplicating filter logic
  #   - Automatic dependency tracking
  #   - Efficient re-computation (only when needed)
  #
  # Filter Logic:
  #   - Start with all variants
  #   - Apply each filter if checkbox is checked
  #   - Filters are cumulative (AND logic)
  #   - Result: Variants matching all selected criteria
  
  filtered_data <- reactive({
    # Start with full dataset
    data <- variant_data
    
    # Filter 1: Biorepository samples
    # Only applies if checkbox is checked
    if (input$biorepository_filter) {
      # Keep only rows where biorepository column has data
      # !is.na() checks for non-missing values
      # != "" checks for non-empty strings
      # Why both? Some missing data is NA, some is empty string
      data <- data %>% filter(!is.na(biorepository) & biorepository != "")
    }
    
    # Filter 2: iPSC cell lines
    if (input$iPSC_line_filter) {
      # Keep only rows where iPSC.line column has data
      # Note: Column name has dot (iPSC.line) not underscore
      data <- data %>% filter(!is.na(iPSC.line) & iPSC.line != "")
    }
    
    # Filter 3: Mouse models
    if (input$mouse_line_filter) {
      # Keep only rows where mouse.line column has data
      data <- data %>% filter(!is.na(mouse.line) & mouse.line != "")
    }
    
    # Return filtered data
    # Result depends on which checkboxes are checked:
    #   - None checked: All 153 variants
    #   - One checked: Subset with that resource
    #   - Multiple checked: Intersection (variants with ALL resources)
    data
  })
  
  # ============================================================
  # DYNAMIC UI GENERATION
  # ============================================================
  
  #' Generate Variant Type Buttons
  #'
  # @description
  # Creates action buttons dynamically based on filtered variant types.
  #
  # Why dynamic?
  #   - Number and types of variants change with filters
  #   - Button list should reflect what's actually displayable
  #   - Automatic adaptation to data changes
  #
  # renderUI() creates UI elements reactively
  #   - Re-runs when filtered_data() changes
  #   - Updates button list in real-time
  
  output$trackButtons <- renderUI({
    
    # Full list of variant types after filtering
    types <- c(
      "missense",
      "nonsense",
      "splice",
      "frameshift",
      "indel",
#remove to show button "cnv",
      "vus",
      "other"
    )
    
    # Exclude problematic types
    # Why exclude?
    #   - "gain exon 3" might be malformed or edge case
    #   - "extra copy" is structural variant without coordinates
    # Better to skip than show broken tracks
    types <- types[types != "gain exon 3"]
    types <- types[types != "extra copy"]
    
    # Create one button per variant type
    # lapply returns a list of button widgets
    #
    # Display label mapping: maps internal variant type names (used as track
    # identifiers in IGV and throughout the data pipeline) to the user-facing
    # button labels shown in the sidebar.
    # The inputId always uses the internal name (via make.names()) so all
    # downstream observers continue to work without modification.
    track_label_map <- c(
      "missense-VUS" = "VUS (missense only)"
    )
    
    buttons <- lapply(types, function(type) {
      # Create action button
      # ID: Unique identifier using variant type
      #     make.names() ensures valid R names (no spaces/special chars)
      #     Example: "missense-VUS" → "missense.VUS"
      # Label: Display text (what user sees)
      #        Uses track_label_map for renamed types, raw type name otherwise
      # as.character() is required: variant.type is a factor, and using a factor
      # as a [[ key makes R treat it as an integer index rather than a string,
      # which causes "subscript out of bounds" for any type not in the map.
      display_label <- if (as.character(type) %in% names(track_label_map)) track_label_map[[as.character(type)]] else as.character(type)
      actionButton(inputId = paste0("addTrack_", make.names(type)), 
                   label = display_label)
    })
    
    # Convert list of buttons to tagList
    # Why tagList? Shiny's way of combining multiple UI elements
    # Renders as a vertical stack of buttons
    do.call(tagList, buttons)
  })
  
  # ============================================================
  # EVENT HANDLING: Button Clicks
  # ============================================================
  
  #' Handle Track Addition Button Clicks
  #'
  # @description
  # Responds to user clicking variant type buttons by adding IGV tracks.
  #
  # Challenge: How to observe dynamically created buttons?
  # Solution: observe() with loop over all variant types
  #
  # Design Pattern: Closure with local()
  #   - Loop variable 'type' would be shared across iterations
  #   - local() creates isolated scope
  #   - Each observeEvent gets its own copy of 't'
  #   - Prevents bugs from variable sharing
  
  observe({
    # Get all variant types (not just filtered ones)
    # Why all? Button IDs exist even when hidden by filters
    types <- c(
      "missense",
      "nonsense",
      "splice",
      "frameshift",
      "indel",
#remove to show button "cnv",
      "vus",
      "other"
    )
    
    # Create observer for each variant type
    for (type in types) {
      # local() creates isolated scope to avoid closure issues
      local({
        # Copy loop variable to local scope
        # Each iteration gets its own 't' variable
        # Without this, all observers would share the same 'type'
        t <- type
        
        # Create observer for this specific button
        # Fires when button is clicked
        observeEvent(input[[paste0("addTrack_", make.names(t))]], {
          
          # Update state: Mark this track as added
          # Get current state
          added_tracks_list <- added_tracks()
          # Set this type to TRUE
          added_tracks_list[[t]] <- TRUE
          # Save updated state
          added_tracks(added_tracks_list)
          
          # Get data for this variant type
          # Uses filtered_data() so respects active filters
          track_data <- filtered_data() %>% filter(variant.type == t)
          
          # Check if we have any data to display
          if (nrow(track_data) > 0) {
            # YES - variants to show
            
            # Convert to GFF3 format
            gff3_data <- create_gff3_data(track_data, variant_type = t)
            
            # Load track into IGV
            # Why loadGFF3TrackFromLocalData?
            #   - Loads from R data frame (not file)
            #   - Fast (no disk I/O)
            #   - Dynamic (updates in real-time)
            loadGFF3TrackFromLocalData(
              session,                         # Shiny session object
              id = "igvShiny_0",               # IGV widget ID
              trackName = t,                   # Track name (e.g., "missense")
              tbl = gff3_data,                 # Data in GFF3 format
              colorTable = color_table,        # Color scheme
              colorByAttribute = "type",       # Color by variant type
              trackHeight = 40,                # Height in pixels
              displayMode = "EXPANDED",        # Show all variants
              visibilityWindow = 1e8           # Show at all zoom levels (100Mb)
            )
          } else {
            # NO - filtered out all variants of this type
            # Load empty track to show "no data" rather than error
            # Why load empty track?
            #   - User feedback (shows track exists but filtered out)
            #   - Prevents confusion (track name visible, just no data)
            #   - Consistent behavior (all clicks produce visible result)
            loadGFF3TrackFromLocalData(
              session,
              id = "igvShiny_0",
              trackName = t,
              tbl = data.frame(),              # Empty data frame
              colorTable = color_table,
              colorByAttribute = "type",
              trackHeight = 40,
              displayMode = "EXPANDED",
              visibilityWindow = 1e8
            )
          }
          
          # Re-center view on SYNGAP1 gene after adding track
          # Why? Ensures user sees the gene region, not random chromosome region
          # "SYNGAP1" is a gene symbol - IGV knows where it is
          showGenomicRegion(session, id = "igvShiny_0", region = "SYNGAP1")
          
          # ── ClinVar pinning ──────────────────────────────────────────────
          # IGV adds new tracks to the bottom of the stack.  If ClinVar is
          # already visible when a variant-type track is introduced for the
          # first time, the new track lands below ClinVar, pushing ClinVar
          # up.  Repinning removes and re-adds ClinVar so it always stays at
          # the very bottom regardless of which variant tracks are enabled.
          if (clinvar_track_added()) {
            repin_clinvar_to_bottom(session, clinvar_data, input$clinvar_size_kb)
          }
          
        }, ignoreInit = TRUE)  # Don't fire on initial app load, only on clicks
      })
    }
  })
  
  # ============================================================
  # CLINVAR TRACK BUTTON HANDLER
  # ============================================================
  
  #' Handle ClinVar Track Button Click
  #'
  # @description
  # Loads the ClinVar variant track into IGV when the user clicks "ClinVar".
  #
  # Design Notes:
  #   - Handled separately from the dynamic variant-type loop because ClinVar
  #     is an external reference dataset, not a filtered subset of internal data.
  #   - Not affected by the biorepository / cell-line / mouse-line checkboxes.
  #   - Uses pre-parsed clinvar_data (built at startup) for instant response.
  #   - Passes the current slider value to create_clinvar_gff3_data() so the
  #     initial load already respects whatever size threshold is set.
  #   - Sets clinvar_track_added(TRUE) so the slider observer knows to react
  #     to future threshold changes.
  #   - Track name "ClinVar" is used so it appears as a labeled track in IGV.
  
  observeEvent(input$addTrack_clinvar, {
    
    # Mark ClinVar track as active so the slider observer can reload it
    clinvar_track_added(TRUE)
    
    # Build GFF3 from parsed ClinVar data, applying the current size threshold.
    # as.numeric() is required because sliderTextInput() returns a character.
    clinvar_gff3 <- create_clinvar_gff3_data(clinvar_data,
                                             max_size_kb = as.numeric(input$clinvar_size_kb))
    
    if (nrow(clinvar_gff3) > 0) {
      # Load ClinVar track into IGV
      # trackName = "ClinVar" sets the label visible in the IGV track panel
      loadGFF3TrackFromLocalData(
        session,
        id                = "igvShiny_0",
        trackName         = "ClinVar",
        tbl               = clinvar_gff3,
        colorTable        = color_table,
        colorByAttribute  = "type",       # "clinvar" type → teal (#00BCD4)
        trackHeight       = 50,           # Slightly taller than internal tracks
        displayMode       = "EXPANDED",   # Show all variants (not collapsed)
        visibilityWindow  = 1e8           # Visible at all zoom levels
      )
    } else {
      # Fallback: load empty track if all rows lacked GRCh38 coordinates
      # Unlikely with a fresh ClinVar download, but handles gracefully
      loadGFF3TrackFromLocalData(
        session,
        id                = "igvShiny_0",
        trackName         = "ClinVar",
        tbl               = data.frame(),
        colorTable        = color_table,
        colorByAttribute  = "type",
        trackHeight       = 50,
        displayMode       = "EXPANDED",
        visibilityWindow  = 1e8
      )
    }
    
    # Re-center view on SYNGAP1 after adding track
    showGenomicRegion(session, id = "igvShiny_0", region = "SYNGAP1")
    
  }, ignoreInit = TRUE)
  
  # ============================================================
  # CLINVAR SIZE SLIDER HANDLER
  # ============================================================
  
  #' Reload ClinVar Track When Size Threshold Changes
  #'
  # @description
  # Responds to changes in the "Max ClinVar variant size" slider by reloading
  # the ClinVar track with only variants at or below the new KB threshold.
  #
  # Only fires if ClinVar track has already been added (clinvar_track_added()
  # is TRUE) — avoids a no-op reload on app startup when the slider is
  # initialized but no track exists yet.
  #
  # Design: Same reload pattern as the filter-change handler used for internal
  # variant tracks — replace the existing track with a freshly filtered one.
  # IGV handles the swap without flicker.
  
  observeEvent(input$clinvar_size_kb, {
    
    # Only reload if the ClinVar track is actually displayed
    if (!clinvar_track_added()) return()
    
    # Re-build GFF3 with the updated size threshold.
    # as.numeric() is required because sliderTextInput() returns a character.
    clinvar_gff3 <- create_clinvar_gff3_data(clinvar_data,
                                             max_size_kb = as.numeric(input$clinvar_size_kb))
    
    if (nrow(clinvar_gff3) > 0) {
      loadGFF3TrackFromLocalData(
        session,
        id                = "igvShiny_0",
        trackName         = "ClinVar",
        tbl               = clinvar_gff3,
        colorTable        = color_table,
        colorByAttribute  = "type",
        trackHeight       = 50,
        displayMode       = "EXPANDED",
        visibilityWindow  = 1e8
      )
    } else {
      # All ClinVar variants exceed the threshold — show empty track
      loadGFF3TrackFromLocalData(
        session,
        id                = "igvShiny_0",
        trackName         = "ClinVar",
        tbl               = data.frame(),
        colorTable        = color_table,
        colorByAttribute  = "type",
        trackHeight       = 50,
        displayMode       = "EXPANDED",
        visibilityWindow  = 1e8
      )
    }
    
  }, ignoreInit = TRUE)  # Don't fire on startup before a track exists
  
  # ============================================================
  # FILTER CHANGE HANDLING
  # ============================================================
  
  #' Update Tracks When Filters Change
  #'
  # @description
  # Refreshes visible tracks when user changes filter checkboxes.
  #
  # Challenge: Filters change data, but tracks are already displayed
  # Solution: Observe filtered_data() and reload all visible tracks
  #
  # Behavior:
  #   - Check biorepository filter → Tracks update to show only those variants
  #   - Uncheck filter → Tracks update to show all variants again
  #
  # Design Decision: Reload vs Hide
  #   - Could hide filtered-out variants (faster)
  #   - Instead, reload entire track (cleaner)
  #   - Trade: Performance for correctness
  #   - Dataset is small (153 variants) so reload is fast enough
  
  observeEvent(filtered_data(), {
    # Get list of currently displayed tracks
    added_tracks_list <- added_tracks()
    
    # Update each displayed track
    # Only updates tracks that user has added (not all possible tracks)
    for (t in names(added_tracks_list)) {
      # Check if this track is supposed to be visible
      if (added_tracks_list[[t]]) {
        # Get filtered data for this variant type
        track_data <- filtered_data() %>% filter(variant.type == t)
        
        # Check if filtered data contains any variants
        if (nrow(track_data) > 0) {
          # YES - show filtered variants
          gff3_data <- create_gff3_data(track_data, variant_type = t)
          
          # Reload track with new data
          # This replaces the existing track
          # IGV handles the update smoothly (no flicker)
          loadGFF3TrackFromLocalData(
            session,
            id = "igvShiny_0",
            trackName = t,
            tbl = gff3_data,
            colorTable = color_table,
            colorByAttribute = "type",
            trackHeight = 40,
            displayMode = "EXPANDED",
            visibilityWindow = 1e8
          )
        } else {
          # NO - all variants filtered out
          # Show empty track
          loadGFF3TrackFromLocalData(
            session,
            id = "igvShiny_0",
            trackName = t,
            tbl = data.frame(),  # Empty
            colorTable = color_table,
            colorByAttribute = "type",
            trackHeight = 40,
            displayMode = "EXPANDED",
            visibilityWindow = 1e8
          )
        }
      }
    }
  }, ignoreInit = TRUE)  # Don't fire on startup, only when filters actually change
  
  # When filters change, all visible variant tracks are reloaded by the
  # observer above. Each reload call appends to the IGV track order, which
  # can push ClinVar up if it was already displayed. Repinning after the
  # reload loop ensures ClinVar always stays at the bottom of the stack.
  observeEvent(filtered_data(), {
    if (clinvar_track_added()) {
      repin_clinvar_to_bottom(session, clinvar_data, input$clinvar_size_kb)
    }
  }, ignoreInit = TRUE)
  
  # ============================================================
  # IGV INITIALIZATION
  # ============================================================
  
  #' Render IGV Genome Browser
  #'
  # @description
  # Initializes the IGV browser widget when app starts.
  #
  # renderIgvShiny() is called once when app loads
  # Creates the IGV instance that will display tracks
  #
  # Configuration:
  #   - Genome: hg38 (GRCh38, current human reference)
  #   - Initial view: SYNGAP1 gene locus
  #   - Display mode: SQUISHED (compact track view)
  #   - Tracks: Empty initially (added by user button clicks)
  
  output$igvShiny_0 <- renderIgvShiny({
    # Parse and validate genome specification
    # genomeName: "hg38" specifies human genome build 38
    # initialLocus: "SYNGAP1" centers view on SYNGAP1 gene
    #   - IGV recognizes gene symbols
    #   - Automatically zooms to gene boundaries
    #   - Alternative: Could use coordinates "chr6:33,420,000-33,448,000"
    genomeOptions <- parseAndValidateGenomeSpec(
      genomeName = "hg38",      # Human genome build 38 (GRCh38)
      initialLocus = "SYNGAP1"  # Start viewing at SYNGAP1 gene
    )
    
    # Create IGV instance
    # displayMode: "SQUISHED" makes tracks compact vertically
    #   - Alternatives: "EXPANDED" (tall), "COLLAPSED" (minimal)
    #   - SQUISHED is good default: Readable but space-efficient
    # tracks: list() starts with no tracks
    #   - User adds tracks by clicking buttons
    #   - Clean initial view (not overwhelming)
    igvShiny(
      genomeOptions,             # Genome configuration
      displayMode = "SQUISHED",  # Compact track display
      tracks = list()            # Start with no tracks
    )
  })
}
# ---------------------------------------------------------------------------
# END SECTION 7
# ---------------------------------------------------------------------------

################################################################################
# SECTION 8: APPLICATION LAUNCH
################################################################################
#
# =============================================================================
# Launch Shiny Application
# =============================================================================
#
# @description
# Starts the web server and opens the application in browser.
#
# shinyApp() combines UI and server into runnable app
# When you run this script, it:
#   1. Loads all the data (CSV reading, API queries, etc.)
#   2. Starts a web server (usually on http://localhost:XXXX)
#   3. Opens your default browser to the app
#   4. Waits for user interaction
#   5. Responds to clicks, filters, etc.
#   6. Stops when you close the browser or press Escape in R console
#
# Run this file with: source("app.R") or click "Run App" in RStudio

shinyApp(ui = ui, server = server)

################################################################################
# END OF APPLICATION
################################################################################
# 
# ==============================================================================
# Summary of Architecture:
# ==============================================================================
#
# 1. DATA PIPELINE
#    Endicott 2026_citizen_variants.csv ─┐
#    2026_citizen_variants.xlsx          ├─ Aggregate counts → Merge → Parse → Coordinate Mapping → Display
#    curesyngap1_variants.csv           ─┘
#
# 2. USER INTERACTIONS
#    Click Button → Filter Data → Format GFF3 → Load Track → Update IGV
#    Check Filter → Filter Data → Reload Tracks → Update IGV
#    Move Size Slider → Re-filter ClinVar by span → Reload ClinVar Track → Update IGV
#
# 3. PERFORMANCE OPTIMIZATIONS
#    - Pre-built caches: svv_build_cache.R fetches and caches both ClinVar
#      data (7-day refresh) and transcript exon structure (6-month refresh)
#    - App reads from cache at startup — no outbound API calls at launch
#    - 2026 Citizen and CureSyngap1 counts aggregated once at startup via table()
#    - Preprocessing: Parse data once at startup, not per interaction
#    - Reactive expressions: Auto-recompute only when inputs change
#
# 4. KEY DESIGN DECISIONS
#    - Separate data processing from UI rendering
#    - Cache building decoupled from app startup (svv_build_cache.R)
#    - Local cDNA→genome arithmetic from cached exon structure
#    - ClinVar data loaded from pre-built cache; refreshed by svv_build_cache.R
#    - Per-database patient counts kept as separate files (not merged) to preserve
#      data provenance and allow independent updates
#    - Use reactive programming for smooth interactions
#    - Generate UI dynamically based on data
#    - Graceful degradation when cache or data file is absent
#    - Clear separation of concerns (functions, sections)
# 5. FUTURE IMPROVEMENTS
#    - Add genome positions directly to CSV (eliminate API dependency)
#    - Add download button for filtered variants
#    - Add variant search functionality
#    - Add zoom to variant feature
#    - Add track reordering
#    - Add patient detail modal on variant click
#
################################################################################
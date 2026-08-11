# SYNGAP1 Variant Viewer

An interactive R Shiny application for visualizing SYNGAP1 genetic variants using the Integrative Genomics Viewer (IGV). This tool enables researchers to explore variant data on a genome browser with dynamic filtering by available research resources.

Developed for **[CURE SYNGAP1](https://curesyngap1.org/)** - a global group of families committed to accelerating the science to cure SYNGAP1 & to supporting each other.

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![R Version](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue)
![Status](https://img.shields.io/badge/status-active-success.svg)

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Architecture](#architecture)
- [File Structure](#file-structure)
- [Performance & Caching](#performance--caching)
- [Maintenance](#maintenance)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Acknowledgments](#acknowledgments)
- [License](#license)

---

## Overview

The SYNGAP1 Variant Viewer is a web-based application designed to help researchers visualize and analyze genetic variants in the SYNGAP1 gene. SYNGAP1 variants are associated with neurodevelopmental disorders, and this tool facilitates:

- **Variant visualization** on chromosome 6 using an embedded genome browser
- **ClinVar integration**: live NCBI ClinVar data fetched automatically and displayed as a toggleable reference track
- **Filtering by research assets** (biorepository samples, iPSC cell lines, mouse models)
- **Interactive exploration** with color-coded variant types
- **Coordinate conversion** from cDNA (transcript-relative) to genomic positions (genome-relative)

The application currently displays **153 variants** with per-variant patient counts drawn from three data sources: the 2024 Citizen Health registry, the 2026 Citizen Health registry, and the CureSyngap1 patient survey.

### Key Technologies

- **R Shiny**: Web application framework
- **igvShiny**: R wrapper for Integrative Genomics Viewer
- **Ensembl REST API**: Transcript exon structure fetch (single call; local arithmetic for all coordinate mapping)
- **NCBI E-utilities**: Live ClinVar variant data retrieval
- **Reactive programming**: Real-time UI updates

---

## Features

### Core Functionality

- ✅ **Interactive Genome Browser**: Embedded IGV displaying variants on chromosome 6
- ✅ **ClinVar Reference Track**: Toggleable track showing all SYNGAP1 ClinVar submissions, fetched live from NCBI and cached for 7 days; always pinned as the lowest track in the browser; variant size filterable via a log-spaced slider (1–50,000 KB, default 100 KB)
- ✅ **Multiple Variant Types**: Missense, nonsense, frameshift, indel, splice VUS, CNV, and other variants
- ✅ **Dynamic Track Loading**: Add/remove variant tracks by type
- ✅ **Color-Coded Display**: Each variant type has a distinct color for easy identification
- ✅ **Smart Filtering**: Filter variants by available research resources
- ✅ **Offline Capability**: Works without internet after initial setup (via caching)
- ✅ **Variant Metadata**: Per-database patient counts (2024 Citizen, 2026 Citizen, CureSyngap1 Survey), research asset availability, and genomic span (Size) displayed in the variant click popup

### Filtering Options

1. **Biorepository samples**: Show only variants with stored biological samples
2. **Cell line availability**: Show only variants with iPSC cell lines
3. **Mouse line availability**: Show only variants with mouse models

Filters are combinatorial (AND logic) - checking multiple filters shows variants that meet ALL criteria.

### Variant Display

Each variant shows:
- Protein change (e.g., "Arg135ter") or cDNA notation (e.g., "c.333del")
- **2024 Citizen Count**: number of patients in the 2024 Citizen Health registry with this variant
- **2026 Citizen Count**: number of patients in the 2026 Citizen Health registry with this variant
- **CureSyngap1 Survey Count**: number of respondents in the CureSyngap1 patient survey with this variant
- Genomic position on chromosome 6
- Variant type classification
- Genomic span (Size), displayed as bp, kb, or Mb depending on variant length

---

## Installation

### Prerequisites

- **R version**: ≥ 4.0.0
- **Operating System**: Windows, macOS, or Linux
- **Internet connection**: Required for initial setup only

### Step 1: Clone the Repository

```bash
git clone https://github.com/svv.intern.dev/SVV_v3_4_0
cd syngap1-variant-viewer
```

### Step 2: Install R Package Dependencies

The app requires several R packages. Install them using the provided script:

```r
# Run the automated installation script
source("svv_packages.R")
```

This will install:
- `shiny` - Web application framework
- `igvShiny` - IGV genome browser wrapper (from Bioconductor)
- `dplyr` - Data manipulation
- `stringr` - String processing
- `httr` - HTTP requests for API
- `jsonlite` - JSON parsing
- `shinyWidgets` - Extended UI widgets (log-spaced ClinVar size slider)
- `readxl` - Read Excel files (2026 Citizen variant registry)

### Step 3: Build the Data Caches

Before launching the app for the first time, run the cache builder to pre-fetch
ClinVar data from NCBI and the transcript exon structure from Ensembl:

```r
source("svv_build_cache.R")
```

This creates:
- `cache/clinvar/clinvar_data.rds` — ~30–60 seconds (NCBI fetch)
- `cache/transcript_structure/ENST00000418600.rds` — ~500 ms (Ensembl fetch)

Cache directories are created automatically by this script. You only need to
re-run it when you want to force a refresh before the automatic expiry windows
(7 days for ClinVar, 6 months for the transcript structure).

### Step 4: Run the Application

```r
# In R console
shiny::runApp("svv_app.R")

# Or in RStudio
# Open svv_app.R and click "Run App" button
```

The app will open in your default web browser. With warm caches the app starts
in under 5 seconds and makes no outbound API calls at launch.

---

## Usage

### Starting the Application

1. Ensure all dependencies are installed (run `svv_packages.R` if needed)
2. Ensure caches are populated (run `svv_build_cache.R` if not already done)
3. Run: `shiny::runApp("svv_app.R")`
3. The app opens in your web browser at `http://localhost:XXXX`

### Basic Workflow

**Step 1: View the genome browser**
- The IGV browser loads centered on the SYNGAP1 gene
- Reference genome: hg38 (GRCh38)
- Chromosome: chr6

**Step 2: Add variant tracks**
- Click buttons in the left sidebar to add variant type tracks
- Example: Click "missense" to display all missense variants; click "VUS (missense only)" to display missense variants of uncertain significance
- Each track appears as a horizontal layer in IGV
- Click **"ClinVar"** (below a separator) to overlay all SYNGAP1 ClinVar submissions as a teal reference track; this track is pinned as the lowest track in the browser and is independent of the research-asset filters
- Use the **"Max ClinVar variant size"** slider beneath the ClinVar button to control which variants appear by their genomic span. The slider is log-spaced (1–50,000 KB) so you get fine-grained control at the low end (1 KB steps for point variants) and coarser steps at the high end (5,000 KB steps for large CNVs). Default is 100 KB. Moving the slider immediately reloads the track if it is already displayed.

**Step 3: Apply filters (optional)**
- The three research asset filter checkboxes appear **above** the ClinVar separator, grouped with the other internal-data controls
- Check "Has biorepository samples" to show only variants with stored samples
- Check "Cell line available" to show only variants with iPSC lines
- Check "Mouse line available" to show only variants with mouse models
- Tracks update automatically when filters change

**Step 4: Explore variants**
- Zoom in/out using IGV controls
- Click on variants to see details
- Tracks show variant names (protein changes or cDNA notation)

### First Run vs. Subsequent Runs

**Before first run** — build the caches with `svv_build_cache.R` (requires internet):
- Takes ~35–65 seconds total
- Queries NCBI E-utilities to fetch all SYNGAP1 ClinVar entries (~30–60 seconds); result cached to `cache/clinvar/` for 7 days
- Makes a single Ensembl API call to fetch the SYNGAP1 transcript exon structure (~500 ms); cached to `cache/transcript_structure/` for 6 months

**Every app launch** (can be offline once caches exist):
- Takes <5 seconds
- Reads both the ClinVar data and transcript exon structure from disk cache
- No outbound API calls at launch

**Refreshing stale caches** — re-run `svv_build_cache.R`:
- ClinVar cache is considered stale after 7 days; re-running the script re-fetches from NCBI (~30–60 seconds)
- Transcript exon structure cache expires after 6 months; the script re-fetches from Ensembl (~500 ms)

---

## Architecture

### Application Structure

The app follows a three-layer architecture:

```
┌─────────────────────────────────────────┐
│ DATA LAYER                              │
│ • Fetch ClinVar data (NCBI, 7-day cache)│
│ • Load 2024_citizen_variants.csv        │
│ • Load 2026_citizen_variants.xlsx       │
│ • Load curesyngap1_variants.csv         │
│ • Aggregate per-variant patient counts  │
│ • Parse variant nomenclature            │
│ • Convert coordinates (API + cache)     │
│ • Format protein changes                │
└─────────────────┬───────────────────────┘
                  ↓
┌─────────────────────────────────────────┐
│ PRESENTATION LAYER (UI)                 │
│ • Sidebar controls                      │
│ • IGV browser widget                    │
│ • Dynamic variant-type buttons          │
│ • Static ClinVar track button           │
└─────────────────┬───────────────────────┘
                  ↓
┌─────────────────────────────────────────┐
│ LOGIC LAYER (Server)                    │
│ • Reactive filtering                    │
│ • Track management                      │
│ • Event handling (incl. ClinVar button) │
└─────────────────────────────────────────┘
```

### Key Components

**1. ClinVar Data Loading**
- Source: NCBI E-utilities API (esearch + esummary)
- Fetch: All SYNGAP1 ClinVar variation IDs, then full metadata in batches of 200
- Cache: `cache/clinvar/clinvar_data.rds` — expires after 7 days, then re-fetched automatically
- Fallback: If NCBI is unreachable, app uses stale cache (or loads empty track as last resort)
- Display: GFF3 track colored teal (#00BCD4); all ClinVar metadata visible on click in IGV

**2. Coordinate Conversion**
- Input: cDNA notation (e.g., "c.333del")
- Extraction: Numeric coordinate (333)
- Approach: Single Ensembl API call fetches full transcript exon structure; all cDNA→genome mappings computed locally via arithmetic
- Output: Genome position (chr6:33,425,796)
- Caching: Exon structure cached to `cache/transcript_structure/` with 6-month expiry

**3. Protein Nomenclature Formatting**
- Converts single-letter amino acid codes to three-letter codes
- Example: "p.R135X" → "Arg135ter"
- Improves readability for biologists

**4. Reactive Filtering**
- Filter changes trigger automatic data recalculation
- Tracks reload with filtered data
- No manual refresh needed
- Note: ClinVar track is unaffected by research-asset filters (it is a reference dataset, not internal data)

**5. GFF3 Track Generation**
- Converts both internal variant data and ClinVar data to GFF3 format
- IGV-compatible genomic feature format
- Internal tracks include: ID, Name, Description, patient count, Size (genomic span in bp or kb)
- ClinVar track includes: all ClinVar metadata fields (germline classification, review status, conditions, etc.) plus Size (genomic span in bp, kb, or Mb)

### Data Flow

```
NCBI E-utilities (esearch + esummary)
  ↓
Cache to cache/clinvar/ (7-day expiry)
  ↓
[On button click: format as GFF3 → ClinVar track in IGV]

2024_citizen_variants.csv
  ↓
Aggregate 2026 Citizen counts (2026_citizen_variants.xlsx → table())
  ↓
Aggregate CureSyngap1 counts (curesyngap1_variants.csv → table())
  ↓
Extract cDNA coordinates
  ↓
Fetch exon structure from Ensembl (1 API call; cached 6 months)
  ↓
Compute genome positions locally via cDNA→exon arithmetic
  ↓
Merge genome positions + source counts
  ↓
User selects filters
  ↓
Filter variant data
  ↓
Generate GFF3 tracks
  ↓
Display in IGV browser
```

For detailed architecture documentation, see **`svv_architecture.Rmd`**.

---

## File Structure

```
syngap1-variant-viewer/
├── app.R                              # Fully annotated application code
├── svv_build_cache.R                  # Standalone script to pre-build both disk caches
├── svv_architecture.Rmd               # Detailed architecture documentation
├── svv_packages.R                     # Automated dependency installation
├── README.md                          # This file
├── 2024_citizen_variants.csv          # 2024 Citizen Health registry (153 variants, 33 columns)
├── 2026_citizen_variants.xlsx         # 2026 Citizen Health registry (per-patient rows)
├── curesyngap1_variants.csv           # CureSyngap1 patient survey (per-patient rows)
├── cache/
│   ├── transcript_structure/      # Ensembl exon structure cache (6-month expiry, auto-generated)
│   │   ├── ENST00000418600.rds    # Exon data frame for SYNGAP1 transcript
│   │   └── ENST00000418600.rds.timestamp  # ISO timestamp of last fetch
│   └── clinvar/                   # ClinVar data cache (7-day expiry, auto-generated)
│       ├── clinvar_data.rds       # Parsed ClinVar data frame
│       └── last_updated.txt       # ISO timestamp of last successful fetch
└── LICENSE                        # License file
```

### File Descriptions

**app.R**
- Main application with extensive inline comments (500+ lines of documentation)
- Production-ready code with explanations of every function and design decision
- Use this file to run the application
- Excellent resource for understanding the codebase

**svv_build_cache.R**
- Standalone script to pre-fetch and cache both external datasets before app launch
- Fetches ClinVar data from NCBI E-utilities and the transcript exon structure from Ensembl
- Writes results to `cache/clinvar/` and `cache/transcript_structure/`
- Re-run to force a refresh before the automatic expiry windows (7 days / 6 months)
- Includes [OK] / [WARN] / [FAIL] diagnostic output, useful for confirming network access and write permissions on a server
- Run once during initial setup; then re-run as needed (or on a scheduled basis on a server)

**svv_architecture.Rmd**
- R Markdown document explaining the overall architecture
- Design philosophy and development timeline
- Key design decisions and tradeoffs
- Common programming patterns used
- Read this to understand the big picture

**svv_packages.R**
- Automated script to install all required R packages
- Handles both CRAN and Bioconductor dependencies
- Run once during initial setup

**2024_citizen_variants.csv**
- 2024 Citizen Health variant registry
- 153 variants across 33 data columns; one row per unique variant
- Includes: variant notation, type, pre-aggregated patient count, research assets
- **Required** for application to function

**2026_citizen_variants.xlsx**
- 2026 Citizen Health variant registry
- Per-patient rows (same variant appears multiple times if multiple patients carry it)
- Columns: cDNA address, protein address, PLP/VUS classification
- Patient counts per variant are derived at startup by counting row occurrences
- **Required** for 2026 Citizen Count display in the variant popup

**curesyngap1_variants.csv**
- CureSyngap1 patient survey export
- Per-patient rows (same variant can appear multiple times)
- Columns: `variant` (cDNA), `protein`
- Patient counts per variant are derived at startup by counting row occurrences
- **Required** for CureSyngap1 Survey Count display in the variant popup

---

## Performance & Caching

The app maintains two separate caches to avoid redundant network calls.

### Ensembl Transcript Structure Cache (6-Month Expiry)

Instead of querying the Ensembl `/map/cdna/` endpoint once per variant coordinate (~43 seconds for 143 variants), the app makes a **single** call to fetch the complete exon structure for the SYNGAP1 transcript, then computes all cDNA→genome mappings locally via arithmetic. The exon structure is cached to disk and auto-refreshes every 6 months to stay current with Ensembl annotation releases.

**How it works:**
1. On first run, `fetch_transcript_exons("ENST00000418600")` retrieves the exon table (~500 ms, 1 API call)
2. Result saved to `cache/transcript_structure/ENST00000418600.rds` with a timestamp
3. `build_cdna_position_map()` converts all coordinates locally (~1 ms) — no further API calls
4. Cache is re-used on subsequent runs until 6 months have elapsed, then auto-refreshed

**To force an immediate transcript cache refresh (e.g., after a known Ensembl release):**
```r
unlink("cache/transcript_structure/ENST00000418600.rds")
unlink("cache/transcript_structure/ENST00000418600.rds.timestamp")
# Restart the app — it will fetch a fresh exon structure
```

### ClinVar Cache (7-Day Expiry)

ClinVar submissions change regularly, so `svv_build_cache.R` re-fetches from NCBI when the cache is older than 7 days.

**How it works:**
1. On startup, the app checks `cache/clinvar/last_updated.txt`
2. If cache is < 7 days old, `clinvar_data.rds` is loaded instantly (<1 second)
3. If cache is absent or expired, the app fetches from NCBI (esearch + esummary, ~30-60 seconds), then saves to `cache/clinvar/`
4. If the fetch fails but a stale cache exists, the stale cache is used with a warning

**To force a re-fetch before 7 days:**
```r
unlink("cache/clinvar/last_updated.txt")
# Restart the app — it will treat the cache as expired
```

### Performance Metrics

| Operation | Old Approach | New Approach |
|-----------|-------------|--------------|
| Ensembl API calls per startup | 143 (one per coordinate) | 1 (exon structure) |
| Ensembl network time (cold) | ~43 seconds | ~500 ms |
| Coordinate mapping time | included above | ~1 ms (local arithmetic) |
| ClinVar fetch (~1,900 variants) | ~30-60 seconds | ~30-60 seconds (unchanged) |
| **Total first-run startup** | **~70–120 seconds** | **~35–65 seconds** |
| Cached startup | <2 seconds | <5 seconds |

**When to clear the ClinVar cache:**
- To force an immediate re-fetch before the 7-day window expires

```r
# Force ClinVar re-fetch on next startup
unlink("cache/clinvar/last_updated.txt")
```

**Cache sizes:**
- Ensembl: ~2KB per coordinate × 143 coordinates ≈ 300KB total
- ClinVar: ~1-3MB for the full SYNGAP1 entry set (single .rds file)

---

## Maintenance

### Updating Variant Data

To update the 2024 Citizen dataset:

1. Export new data from your master database (Excel)
2. Ensure column names match the current structure
3. Replace `2024_citizen_variants.csv`
4. Restart the application

To update the 2026 Citizen dataset:

1. Replace `2026_citizen_variants.xlsx` with the new export
2. Restart the application — patient counts are re-aggregated automatically at startup

To update the CureSyngap1 survey data:

1. Replace `curesyngap1_variants.csv` with the new export
2. Restart the application — patient counts are re-aggregated automatically at startup

The app will automatically:
- Parse new variants
- Map any new cDNA coordinates using the cached exon structure (no additional API calls)
- Re-aggregate patient counts from the updated source files
- Update the visualization

### Updating ClinVar Data

ClinVar data is refreshed automatically every 7 days on startup — no manual action is needed. To force an immediate refresh:

```r
unlink("cache/clinvar/last_updated.txt")
# Then restart the app
```

### Adding New Variant Types

If your data includes new variant classifications:

1. Add to the variant type normalization `case_when` block in Section 6 of `app.R`
2. Add a color entry to `color_table` (also in Section 6, just below the normalization block)
3. No other changes needed — app generates buttons dynamically

Example:
```r
# In variant type normalization
variant.type = case_when(
  variant.type %in% c('splice site') ~ 'splice',
  # ... existing mappings ...
)

# In color table
color_table <- list(
  splice = "#FF5733",  # New color for splice variants
  # ... existing colors ...
)
```

### Changing Genome Build

Currently configured for **hg38**. To switch to hg19:

1. Update `genomeName = "hg19"` in the `renderIgvShiny` block (~line 1385 in `svv_app.R`)
2. Clear transcript cache: `unlink("cache/transcript_structure/", recursive = TRUE)`
3. Restart app (will re-fetch exon structure and remap coordinates for the new build)

Note: The ClinVar GRCh38 coordinate filter in `create_clinvar_gff3_data` would also need updating if switching builds.

### Changing Transcript

Currently uses **ENST00000418600** (SYNGAP1 canonical transcript).

To change:
1. Update the transcript ID in `fetch_transcript_exons()` (in `svv_build_cache.R`) and the `TRANSCRIPT_CACHE_FILE` constant (~line 449 in `svv_app.R`)
2. Clear transcript cache: `unlink("cache/transcript_structure/", recursive = TRUE)`
3. Restart app — it will fetch the exon structure for the new transcript and remap all coordinates

---

## Troubleshooting

### Common Issues

**Problem: `Error in library(igvShiny) : there is no package called 'igvShiny'`**

**Solution:**
```r
source("svv_packages.R")
# Or manually:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("igvShiny")
```

---

**Problem: ClinVar track is empty or missing**

**Solution:** This can happen if the NCBI fetch failed on startup and no cache existed. Check the R console for ClinVar-related warnings. Try:
```r
unlink("cache/clinvar/last_updated.txt")
# Then restart the app
```
If the problem persists, verify your internet connection and that `https://eutils.ncbi.nlm.nih.gov` is reachable.

---

**Problem: ClinVar track shows outdated data**

**Solution:** The cache refreshes automatically every 7 days. To force an immediate re-fetch:
```r
unlink("cache/clinvar/last_updated.txt")
# Restart the app
```

---

**Problem: App loads slowly or takes ~35–65 seconds on first launch**

**Solution:** The slow startup is caused by cache-building, which should be run
separately via `svv_build_cache.R` *before* launching the app:

```r
source("svv_build_cache.R")
```

Once both caches are populated, the app starts in under 5 seconds and makes no
outbound API calls at launch.

---

**Problem: No variants displayed on tracks**

**Solution:** Check if filters are too restrictive. Uncheck all filter boxes to see all variants, then reapply filters incrementally.

---

**Problem: Variants missing genome positions**

**Solution:** Some variants (structural, intronic) may not have valid cDNA coordinates. These show as NA in positions and are skipped in visualization. This is expected behavior.

---

**Problem: Cache directory not found**

**Solution:**
```r
dir.create("cache/transcript_structure", recursive = TRUE)
dir.create("cache/clinvar", recursive = TRUE)
```

---

### Logging & Debugging

Enable console logging to track API calls:

**`svv_build_cache.R` console output** (run before first launch):
- `"ClinVar: querying NCBI esearch for SYNGAP1 variation IDs..."` — ClinVar fetch started
- `"ClinVar: found X variation IDs."` — esearch succeeded
- `"ClinVar: fetching summaries in N batches..."` — esummary batches in progress
- `"ClinVar: successfully parsed X variant records."` — ClinVar data ready
- `"[OK]   ClinVar data loaded: X variants"` — cache written successfully
- `"Ensembl: fetching exon structure for ENST00000418600..."` — transcript fetch started
- `"Ensembl: retrieved X exons (strand -1)"` — exon structure ready
- `"[OK]   Transcript exon structure fetched and cached: X exons"` — cache written successfully

**`svv_app.R` console output** (at each app launch):
- `"ClinVar: loading from cache..."` — cache present and loaded
- `"Loading transcript exon structure from cache..."` — cache present and loaded
- `"cDNA coordinate X is out of transcript range..."` — coordinate exceeds transcript length (expected for some structural variants)

**`svv_app.R` warnings to watch for** (indicate missing cache — run `svv_build_cache.R`):
- `"ClinVar cache not found..."` — ClinVar track will be empty this session
- `"Transcript cache not found..."` — all genome positions will be NA this session

These messages are informational and don't indicate problems unless all variants fail.

---

## Citation

### Software Citations

If you use this tool in your research, please cite:

**This application:**
```
SYNGAP1 Variant Viewer
Developed for Cure SYNGAP1
https://github.com/cmcneil-02/syngap1-variant-viewer
```

**igvShiny (required dependency):**
```
Shannon P, Gladki A, Scigocka K (2024). igvShiny: igvShiny: a wrapper of 
Integrative Genomics Viewer (IGV - an interactive tool for visualization 
and exploration integrated genomic data). R package version 1.2.0, 
https://gladkia.github.io/igvShiny/, https://github.com/gladkia/igvShiny
```

**Ensembl REST API:**
```
Yates AD, et al. (2020). Ensembl 2020. Nucleic Acids Research, 
48(D1), D682-D688. doi: 10.1093/nar/gkz966
```

**NCBI E-utilities (ClinVar data):**
```
Sayers EW, et al. (2022). Database resources of the National Center for 
Biotechnology Information. Nucleic Acids Research, 50(D1), D20-D26. 
doi: 10.1093/nar/gkab1112
```

### BibTeX Entries

```bibtex
@software{syngap1_viewer,
  title = {SYNGAP1 Variant Viewer},
  author = {CURE SYNGAP1 Team},
  organization = {Cure SYNGAP1},
  year = {2025},
  url = {https://github.com/cmcneil-02/syngap1-variant-viewer}
}

@Manual{igvShiny,
  title = {igvShiny: igvShiny wrapper of Integrative Genomics Viewer},
  author = {Paul Shannon and Arkadiusz Gladki and Karolina Scigocka},
  year = {2024},
  note = {R package version 1.2.0},
  url = {https://github.com/gladkia/igvShiny}
}
```

---

## Acknowledgments

### Organizations

**Cure SYNGAP1**
- Website: https://curesyngap1.org/
- Mission: Finding a cure for SYNGAP1-related disorders
- This application was developed to support SYNGAP1 research and patient advocacy efforts

**Citizen Health**
- Patient variant data source and registry partner

### Software & Tools

- **igvShiny developers**: Paul Shannon, Arkadiusz Gladki, Karolina Scigocka
- **Ensembl**: Genome annotation and REST API services
- **NCBI**: ClinVar database and E-utilities API
- **R Shiny team**: Application framework
- **SYNGAP1 research community**: Collaborative data sharing and research efforts

### Contributors

- Collin McNeil
- Chloe Kaufman

---

## License

- MIT License - See LICENSE file for details

This software is provided for research and academic use.

---

## Contact & Support

**Issues & Questions:**
- Open an issue on GitHub: https://github.com/cmcneil-02/SYNGAP1-Variant-Viewer/issues

**CURE SYNGAP1:**
- Website: https://curesyngap1.org/

**Collaboration Inquiries:**
- [Contact method to be determined]

---

## Data Source

Variant data is drawn from three sources:

**Citizen Health (2024 registry)** — `2024_citizen_variants.csv`
- 153 unique SYNGAP1 variants, one row per variant
- Pre-aggregated patient count per variant
- Research asset availability (biorepository samples, cell lines, mouse models)
- Clinical and functional annotations

**Citizen Health (2026 registry)** — `2026_citizen_variants.xlsx`
- Per-patient rows; patient count per variant derived at startup
- 263 total entries across 196 unique variants

**CureSyngap1 Patient Survey** — `curesyngap1_variants.csv`
- Per-patient rows; patient count per variant derived at startup
- 287 total entries across 202 unique variants

All three sources are loaded at app startup. Per-variant counts from the 2026 Citizen and CureSyngap1 datasets are displayed alongside the 2024 Citizen count in the variant click popup.

Data is publicly available through Citizen Health and Cure SYNGAP1.

---

## Version History

## v3.4.0 (2026‑08‑14)
- Overhauled internal variant color system: Replaced legacy categories (missense-VUS, gain, loss, intronic) with a new, simplified set: missense, nonsense, splice, cnv, other. Updated global color table, IGVShiny track attributes, GFF3 formatting, and all server/UI references. Removed unused historical categories and ensured each variant type maps cleanly to a single color
- Rebuilt internal GFF3 generator for consistency: Standardized attribute fields (Name, VariantType, Size, SourceCounts, Color). Cleaned attribute ordering and removed deprecated fields. Ensured all internal variants produce valid IGV‑compatible GFF3 entries with correct color attributes
- Fixed coordinate mismatches between cDNA and genomic positions: Normalized HGVS identifiers and cleaned variant nomenclature. Corrected several inherited cDNA→genome arithmetic inconsistencies. Ensured all internal variants align correctly with transcript structure.
- Improved IGVShiny track reliability: Rewrote portions of track‑formatting logic to prevent silent failures caused by malformed attributes. Ensured internal tracks load with correct colors, names, and popup metadata when IGVShiny renders normally.
- Major Shiny architecture cleanup and modernization: Reorganized app files into a clearer, more maintainable structure. Added detailed comments, consistent formatting, and improved readability across global.R, server.R, and helper modules. Removed dead code, unused variables, and outdated logic paths
- General stability improvements: Fixed small bugs affecting variant filtering, popup display, and reactive track updates. Improved error‑handling around data loading and preprocessing. No changes to ClinVar caching or external API behavior

### v3.3.0 (2026-04-05)
- Added **per-database patient count display** in the variant click popup — each variant now shows three separate count fields: **2024 Citizen Count**, **2026 Citizen Count**, and **CureSyngap1 Survey Count**, replacing the previous single "Number of patients" field
- Added two new data sources: `2026_citizen_variants.xlsx` (2026 Citizen Health registry) and `curesyngap1_variants.csv` (CureSyngap1 patient survey); both are per-patient-row files whose counts are aggregated via `table()` at startup
- Added `library(readxl)` dependency to read the 2026 Citizen `.xlsx` file; `svv_packages.R` updated to install and verify `readxl` automatically
- Both new loaders wrapped in `tryCatch` — if either file is missing the app starts normally with `0` for that source's count rather than crashing
- Renamed data files for clarity: `updatedCitizen191.csv` → `2024_citizen_variants.csv`, `curesyngap1_variants_clean.csv` → `curesyngap1_variants.csv`; all code and documentation updated accordingly
- No changes to coordinate mapping, ClinVar integration, filters, UI layout, or caches; no cache invalidation required

### v3.2.1 (2026-04-03)
- Fixed **Size** display for single-base variants (insertions, deletions, SNVs, etc.) showing **"0 bp"** instead of **"1 bp"** — affected both internal variants and ClinVar variants
- Root cause: both `create_gff3_data()` and `create_clinvar_gff3_data()` computed size as `end - start` (or `clinvar_end - clinvar_start`), which equals `0` when start and end coordinates are identical (point variants); the v3.0.1 fix that removed the ClinVar `+1` offset meant ClinVar SNVs were also susceptible
- Fix: Replaced all raw difference expressions with `pmax(..., 1L)` so the size floors at 1 bp; applied across both branches of the internal attributes `ifelse` and all three tiers of the ClinVar nested `ifelse` (bp / kb / Mb)
- No UI, filter, data pipeline, or cache changes; no cache invalidation required

### v3.2.0 (2026-03-23)
- Added **Size** field to the feature popup for all variant tracks (internal and ClinVar)
- Internal variants display size as bp (< 1,000 bp) or kb (≥ 1,000 bp), computed from genomic start/end coordinates
- ClinVar variants display size as bp, kb, or Mb depending on span; three-tier format accounts for the wide range of ClinVar entries (SNVs through multi-megabase CNVs)
- No UI, filter, data pipeline, or API changes; no cache invalidation required

### v3.1.0 (2026)
- Separated cache-building code from the main app into `svv_build_cache.R`; the app now loads both caches from disk at startup rather than fetching from NCBI/Ensembl on each launch
- Reorganized sidebar: research asset filter checkboxes (biorepository, cell line, mouse line) moved above the horizontal separator so all internal-data controls are grouped together, with the ClinVar section clearly separated below
- Renamed "missense-VUS" track button label to "VUS (missense only)" for clarity; internal track identifier and all data pipeline references are unchanged
- ClinVar track is now pinned as the lowest track in the IGV browser — enabling ClinVar first and then adding a variant-type track, or changing filters while ClinVar is active, will no longer push ClinVar above newly loaded tracks

### v3.0.1 (2026)
- Fixed bug where ClinVar single-base variants (SNVs, single-base deletions) were displayed as 2-base spans in IGV
- Root cause: `parse_esummary_record` was applying a `+1` offset to `end` coordinates whenever `start == stop`, under the incorrect assumption that IGV could not render a zero-width feature; IGV handles `start == end` correctly (as confirmed by the Citizen Health track, which has always used point coordinates without adjustment)
- Fix: Removed the `+1` offset; coordinates are now passed through from NCBI esummary without modification
- No functional, UI, or API changes; no cache invalidation required beyond clearing the ClinVar cache to re-fetch clean coordinates (see Troubleshooting)

### v3.0.0 (2026)
- Replaced per-coordinate Ensembl `/map/cdna/` memoization with a single exon-structure fetch and local cDNA→genome arithmetic
- Cold-cache startup time reduced from ~8–12 minutes to under 30 seconds
- Transcript exon structure cached to `cache/transcript_structure/` with 6-month auto-refresh
- Removed `memoise` and `cachem` dependencies
- Prior approach preserved on the `slow-cache-option` branch for reference

### v2.1.0 (2026)
- Added live ClinVar integration via NCBI E-utilities (esearch + esummary)
- ClinVar data fetched automatically on startup and cached locally for 7 days
- ClinVar variants displayed as a toggleable teal reference track in IGV
- All ClinVar metadata (germline classification, review status, conditions, etc.) accessible on variant click
- Log-spaced size filter slider (1–50,000 KB, default 100 KB) controls which ClinVar variants are shown by genomic span; reloads track live on change

### v1.0.0 (2025)
- Initial release
- Support for 153 variants
- Three filter types (biorepository, iPSC, mouse models)
- IGV browser integration
- Ensembl API with permanent caching
- Color-coded variant types
- Reactive filtering

---

**Community Contributions Welcome!**

---

**Last Updated:** 08-14-2026  
**Repository:** SVV_v3_4_0  
**Documentation Version:** 3.4.0  
**Maintained by:** CURE SYNGAP1 Team

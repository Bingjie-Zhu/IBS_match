# IBSmatch

An R package for IBS (Identity-by-State) based DNA-RNA sample matching.

## Overview

**IBSmatch** provides a complete workflow for verifying that DNA and RNA
samples collected from the same individual have been correctly labelled and
paired.  It takes an IBS similarity matrix (e.g. computed by PLINK's
`--distance ibs` command) and a table of expected DNA-RNA pairings, then:

1. Calculates a Top-1% IBS threshold from the background distribution.
2. Plots the background distribution against the known matched-pair
   distribution.
3. Identifies Reciprocal Best Hit (RBH) matches.
4. (Optionally) audits RBH results against a gold-standard pairing list.

## Installation

```r
# install.packages("remotes")
remotes::install_github("Bingjie-Zhu/IBS_match")
```

## Quick Start

```r
library(IBSmatch)

# 1. Build the similarity sub-matrix from PLINK output
sub_mat <- build_ibs_submatrix(
  dist_path    = "rna_dist.mdist",
  id_path      = "rna_dist.mdist.id",
  sample_pairs = "matched_samples.txt"   # columns: DNA_name, RNA_name
)

# 2. Calculate the top-1% IBS threshold
threshold <- calc_ibs_threshold(sub_mat)
message("Top 1% threshold: ", round(threshold, 4))

# 3. Plot background vs matched-pair IBS distributions
plot_ibs_distribution(sub_mat, title = "Wheat DNA-RNA IBS Distribution")

# 4. Find Reciprocal Best Hit matches (all pairs, no threshold)
rbh <- find_rbh_matches(sub_mat)

# 4b. Optionally filter by IBS score
rbh_filtered <- find_rbh_matches(sub_mat, ibs_threshold = threshold)

# 5. Audit against expected pairings
expected <- data.frame(
  DNA          = rownames(sub_mat),
  Expected_RNA = colnames(sub_mat),
  stringsAsFactors = FALSE
)
audit <- audit_matches(rbh, expected)
write.csv(audit, "audit_result.csv", quote = FALSE, row.names = FALSE)
```

## Exported Functions

| Function | Description |
|---|---|
| `build_ibs_submatrix()` | Read `.mdist` files and build the DNA × RNA similarity sub-matrix |
| `calc_ibs_threshold()` | Compute the 99th-percentile (top-1%) IBS threshold |
| `plot_ibs_distribution()` | ggplot2 density plot of background vs matched-pair IBS scores |
| `find_rbh_matches()` | Reciprocal Best Hit matching with optional IBS threshold filtering |
| `audit_matches()` | Compare RBH results against a gold-standard pairing list |

## Dependencies

- [ggplot2](https://ggplot2.tidyverse.org/)

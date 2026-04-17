#' Find Reciprocal Best Hit (RBH) Matches
#'
#' Identifies DNA-RNA sample pairs that are each other's best IBS match
#' (Reciprocal Best Hit). In the forward direction every DNA sample is paired
#' with the RNA sample that has the highest IBS score. In the reverse direction
#' every RNA sample is paired with the DNA sample that has the highest IBS
#' score. Only pairs that appear in *both* directions are retained.
#'
#' By default **all** RBH pairs are returned. Set `ibs_threshold` to a value
#' between 0 and 1 to keep only pairs whose IBS score meets or exceeds that
#' threshold.
#'
#' @param sub_mat A numeric matrix with DNA samples as rows and RNA samples as
#'   columns. Values should be IBS similarity scores (0-1). Row names must
#'   identify the DNA samples and column names must identify the RNA samples.
#'   `NA` values are treated as non-matches (i.e. the corresponding sample is
#'   never chosen as a best hit).
#' @param ibs_threshold Either `NULL` (default, no filtering) or a single
#'   numeric value in [0, 1]. When provided, only RBH pairs with
#'   `IBS_Score >= ibs_threshold` are returned.
#'
#' @return A data frame with three columns:
#' \describe{
#'   \item{DNA}{Character. DNA sample identifier.}
#'   \item{RNA}{Character. RNA sample identifier.}
#'   \item{IBS_Score}{Numeric. IBS similarity score for the pair.}
#' }
#'
#' @examples
#' \dontrun{
#' sim_mat <- matrix(runif(100), nrow = 10,
#'                   dimnames = list(paste0("DNA", 1:10), paste0("RNA", 1:10)))
#' rbh <- find_rbh_matches(sim_mat)
#' rbh_filtered <- find_rbh_matches(sim_mat, ibs_threshold = 0.8)
#' }
#'
#' @export
find_rbh_matches <- function(sub_mat, ibs_threshold = NULL) {
  if (!is.matrix(sub_mat) || !is.numeric(sub_mat)) {
    stop("'sub_mat' must be a numeric matrix.")
  }
  if (is.null(rownames(sub_mat)) || is.null(colnames(sub_mat))) {
    stop("'sub_mat' must have row names (DNA) and column names (RNA).")
  }
  if (!is.null(ibs_threshold)) {
    if (!is.numeric(ibs_threshold) || length(ibs_threshold) != 1 ||
        ibs_threshold < 0 || ibs_threshold > 1) {
      stop("'ibs_threshold' must be a single numeric value between 0 and 1.")
    }
  }

  # Replace NAs with -1 so they are never chosen as the best hit
  mat <- sub_mat
  mat[is.na(mat)] <- -1

  # Forward: for each DNA row find the RNA column with the highest IBS score
  best_rna_idx <- max.col(mat, ties.method = "first")
  dna_to_rna <- data.frame(
    DNA       = rownames(mat),
    RNA       = colnames(mat)[best_rna_idx],
    IBS_Score = mat[cbind(seq_len(nrow(mat)), best_rna_idx)],
    stringsAsFactors = FALSE
  )

  # Reverse: for each RNA column find the DNA row with the highest IBS score
  t_mat <- t(mat)
  best_dna_idx <- max.col(t_mat, ties.method = "first")
  rna_to_dna <- data.frame(
    RNA = rownames(t_mat),
    DNA = colnames(t_mat)[best_dna_idx],
    stringsAsFactors = FALSE
  )

  # Intersection: keep only pairs that are mutually best hits
  rbh <- merge(dna_to_rna, rna_to_dna, by = c("DNA", "RNA"))

  # Optional threshold filtering
  if (!is.null(ibs_threshold)) {
    rbh <- rbh[rbh$IBS_Score >= ibs_threshold, , drop = FALSE]
  }

  rownames(rbh) <- NULL
  return(rbh)
}

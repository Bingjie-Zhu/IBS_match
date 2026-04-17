#' Audit RBH Matches Against Expected Sample Pairs
#'
#' Compares a set of Reciprocal Best Hit (RBH) matches against a gold-standard
#' list of expected DNA-RNA pairings and annotates each result row with whether
#' the predicted RNA matches the expected RNA.
#'
#' @param rbh_result A data frame returned by \code{\link{find_rbh_matches}},
#'   containing at minimum the columns `DNA`, `RNA`, and `IBS_Score`.
#' @param expected_pairs A data frame with two columns:
#' \describe{
#'   \item{DNA}{Character. DNA sample identifier (must match identifiers used
#'     in `rbh_result$DNA`).}
#'   \item{Expected_RNA}{Character. The RNA sample that is expected to match
#'     each DNA sample according to the gold standard.}
#' }
#'
#' @return A data frame that combines `rbh_result` with `Expected_RNA` from
#'   `expected_pairs` (matched on `DNA`) and an additional logical column
#'   `Is_Match` that is `TRUE` when the predicted `RNA` equals `Expected_RNA`
#'   and `FALSE` otherwise. Rows are sorted so that mismatches appear first
#'   (ascending `Is_Match`), with ties broken by descending `IBS_Score`.
#'
#' @examples
#' \dontrun{
#' rbh <- find_rbh_matches(sim_mat)
#' expected <- data.frame(
#'   DNA          = paste0("DNA", 1:10),
#'   Expected_RNA = paste0("RNA", 1:10),
#'   stringsAsFactors = FALSE
#' )
#' audit <- audit_matches(rbh, expected)
#' }
#'
#' @export
audit_matches <- function(rbh_result, expected_pairs) {
  if (!is.data.frame(rbh_result)) {
    stop("'rbh_result' must be a data frame (e.g. from find_rbh_matches()).")
  }
  if (!all(c("DNA", "RNA", "IBS_Score") %in% names(rbh_result))) {
    stop("'rbh_result' must contain columns: DNA, RNA, IBS_Score.")
  }
  if (!is.data.frame(expected_pairs)) {
    stop("'expected_pairs' must be a data frame.")
  }
  if (!all(c("DNA", "Expected_RNA") %in% names(expected_pairs))) {
    stop("'expected_pairs' must contain columns: DNA, Expected_RNA.")
  }

  audit_result <- merge(rbh_result, expected_pairs, by = "DNA", all.x = TRUE)
  audit_result$Is_Match <- audit_result$RNA == audit_result$Expected_RNA

  audit_result <- audit_result[
    order(audit_result$Is_Match, -audit_result$IBS_Score), ,
    drop = FALSE
  ]

  rownames(audit_result) <- NULL
  return(audit_result)
}

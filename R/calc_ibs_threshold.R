#' Calculate the Top 1% IBS Threshold
#'
#' Computes the 99th-percentile (top 1%) of all pairwise IBS similarity values
#' in a sub-matrix whose rows represent DNA samples and whose columns represent
#' RNA samples.
#'
#' @param sub_mat A numeric matrix with DNA samples as rows and RNA samples as
#'   columns. Values should be IBS similarity scores (0–1).
#'
#' @return A single numeric value: the IBS similarity score at the 99th
#'   percentile of the background distribution.
#'
#' @examples
#' \dontrun{
#' sim_mat <- matrix(runif(100), nrow = 10,
#'                   dimnames = list(paste0("DNA", 1:10), paste0("RNA", 1:10)))
#' threshold <- calc_ibs_threshold(sim_mat)
#' }
#'
#' @export
calc_ibs_threshold <- function(sub_mat) {
  if (!is.matrix(sub_mat) || !is.numeric(sub_mat)) {
    stop("'sub_mat' must be a numeric matrix.")
  }

  all_values <- as.vector(sub_mat)
  threshold <- quantile(all_values, 0.99, na.rm = TRUE)
  return(threshold)
}

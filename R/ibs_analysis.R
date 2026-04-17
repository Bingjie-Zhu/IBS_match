#' Build an IBS Similarity Sub-Matrix from a Distance Matrix File
#'
#' A convenience helper that reads a PLINK-style `.mdist` distance matrix and
#' its accompanying `.mdist.id` sample-ID file, converts distances to
#' similarities (`1 - distance`), and extracts the sub-matrix whose rows
#' correspond to DNA samples and whose columns correspond to RNA samples.
#'
#' @param dist_path Character. Path to the `.mdist` distance matrix file
#'   (whitespace-separated, no header).
#' @param id_path Character. Path to the `.mdist.id` sample-ID file. The
#'   function reads the second column as sample names.
#' @param sample_pairs A data frame (or file path to a whitespace-separated
#'   table) with two columns: DNA sample names (first) and RNA sample names
#'   (second). Column names are ignored; positional order is used.
#'
#' @return A numeric matrix with DNA samples as rows and RNA samples as columns,
#'   containing IBS similarity scores (0-1).
#'
#' @examples
#' \dontrun{
#' sub_mat <- build_ibs_submatrix(
#'   dist_path    = "rna_dist.mdist",
#'   id_path      = "rna_dist.mdist.id",
#'   sample_pairs = "matched_samples.txt"
#' )
#' }
#'
#' @export
build_ibs_submatrix <- function(dist_path, id_path, sample_pairs) {
  dist_mat <- as.matrix(read.table(dist_path))
  ids      <- read.table(id_path)[, 2]

  sim_mat <- 1 - dist_mat
  rownames(sim_mat) <- ids
  colnames(sim_mat) <- ids

  if (is.character(sample_pairs)) {
    sample_pairs <- read.table(sample_pairs)
    }
  dna_names <- as.character(sample_pairs[, 2])
  rna_names <- as.character(sample_pairs[, 3])
  sub_mat <- sim_mat[dna_names, rna_names]
  return(sub_mat)
}

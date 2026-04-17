#' Plot IBS Background and Matched-Pair Distributions
#'
#' Draws a density plot that overlays the background distribution of all
#' pairwise IBS values with the distribution of the known matched pairs
#' (diagonal of `sub_mat`). A vertical dashed line marks the top-1% threshold
#' of the background distribution.
#'
#' @param sub_mat A numeric matrix with DNA samples as rows and RNA samples as
#'   columns. Values should be IBS similarity scores (0–1). Row and column
#'   names must be set and correspond to DNA and RNA sample identifiers,
#'   respectively.
#' @param title Character string used as the plot title.
#'   Defaults to `"DNA-RNA IBS Distribution"`.
#'
#' @return A \code{ggplot} object (invisibly). The plot is also printed to the
#'   active graphics device.
#'
#' @examples
#' \dontrun{
#' sim_mat <- matrix(runif(100), nrow = 10,
#'                   dimnames = list(paste0("DNA", 1:10), paste0("RNA", 1:10)))
#' plot_ibs_distribution(sim_mat)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_density geom_vline annotate labs
#'   scale_fill_manual theme_minimal
#' @export
plot_ibs_distribution <- function(sub_mat,
                                  title = "DNA-RNA IBS Distribution") {
  if (!is.matrix(sub_mat) || !is.numeric(sub_mat)) {
    stop("'sub_mat' must be a numeric matrix.")
  }

  n_diag <- min(nrow(sub_mat), ncol(sub_mat))

  all_pairwise_ibs <- as.vector(sub_mat)
  matched_ibs <- diag(sub_mat[seq_len(n_diag), seq_len(n_diag), drop = FALSE])

  bg_threshold <- quantile(all_pairwise_ibs, 0.99, na.rm = TRUE)

  plot_data <- rbind(
    data.frame(IBS = all_pairwise_ibs, Group = "Background",
               stringsAsFactors = FALSE),
    data.frame(IBS = matched_ibs,      Group = "Matched Pairs",
               stringsAsFactors = FALSE)
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = IBS, fill = Group)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::geom_vline(xintercept = bg_threshold,
                        linetype = "dashed", color = "blue") +
    ggplot2::annotate("text",
                      x = bg_threshold,
                      y = Inf,
                      label = paste0("BG Top 1%\n(",
                                     round(bg_threshold, 4), ")"),
                      hjust = -0.05, vjust = -0.5, color = "blue", size = 3) +
    ggplot2::scale_fill_manual(
      values = c("Background" = "grey", "Matched Pairs" = "#E41A1C")
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = paste("Samples:", n_diag,
                       "| Background Top 1% threshold:",
                       round(bg_threshold, 4)),
      x = "IBS Similarity",
      y = "Density"
    ) +
    ggplot2::theme_minimal()

  print(p)
  invisible(p)
}

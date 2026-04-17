# Function to find Reciprocal Best Hits (RBH) from a similarity matrix

find_rbh_matches <- function(similarity_matrix, ibs_threshold = 0.7) {
  # Check if the similarity matrix is valid
  if (is.null(similarity_matrix) || !is.matrix(similarity_matrix)) {
    stop('Input must be a non-null similarity matrix.')
  }

  # Get the number of rows and columns in the similarity matrix
  n <- nrow(similarity_matrix)
  m <- ncol(similarity_matrix)

  # Ensure that the dimensions are consistent for RBH calculation
  if (n != m) {
    stop('Similarity matrix must be square (same number of rows and columns).')
  }

  # Initialize an empty list for RBH matches
  rbh_matches <- list()

  # Iterate through each element in the similarity matrix
  for (i in 1:n) {
    for (j in 1:n) {
      # Check if the value meets the IBS threshold
      if (similarity_matrix[i, j] >= ibs_threshold) {
        # Check for reciprocal best hit
        if (similarity_matrix[j, i] >= ibs_threshold) {
          rbh_matches[[paste(i, j, sep = '-')]] <- list(
            hit1 = i,
            hit2 = j,
            score1 = similarity_matrix[i, j],
            score2 = similarity_matrix[j, i]
          )
        }
      }
    }
  }

  return(rbh_matches)
}

# Example Usage:
# similarity_matrix <- matrix(c(1, 0.8, 0.4, 0.8, 1, 0.3, 0.2, 0.3, 1), nrow = 3)
# rbh_results <- find_rbh_matches(similarity_matrix, ibs_threshold = 0.7)
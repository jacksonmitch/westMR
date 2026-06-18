row_logsumexp <- function(log_mat) {
  row_max <- apply(log_mat, 1, max)
  row_max + log(rowSums(exp(log_mat - row_max)))
}

# Numerically stable softmax for log weights
row_softmax <- function(log_mat) {
  row_max <- apply(log_mat, 1, max)
  w <- exp(log_mat - row_max)
  w / rowSums(w)
}

common_eta <- function(B, beta = numeric(0)) {
  B <- as.matrix(B)
  beta <- as.numeric(beta)

  if (ncol(B) == 0) {
    return(rep(0, nrow(B)))
  }

  if (length(beta) != ncol(B)) {
    stop(
      "Length of beta must equal ncol(B). Got length(beta) = ",
      length(beta),
      " and ncol(B) = ",
      ncol(B),
      "."
    )
  }

  drop(B %*% beta)
}

linear_predictor_matrix <- function(A, B, beta_g, beta) {
  A <- as.matrix(A)

  if (is.null(B)) {
    B <- matrix(numeric(0), nrow = nrow(A), ncol = 0)
  } else {
    B <- as.matrix(B)
  }

  beta_g <- as.matrix(beta_g)
  beta <- as.numeric(beta)

  sweep(A %*% t(beta_g), 1, common_eta(B, beta), "+")
}

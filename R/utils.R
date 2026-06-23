row_logsumexp <- function(log_mat) {
  log_mat <- as.matrix(log_mat)
  row_max <- row_max_base(log_mat)
  row_max + log(rowSums(exp(log_mat - row_max)))
}

row_softmax <- function(log_mat) {
  log_mat <- as.matrix(log_mat)
  row_max <- row_max_base(log_mat)
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

  eta <- sweep(A %*% t(beta_g), 1, common_eta(B, beta), "+")
  eta <- as.matrix(eta)
  storage.mode(eta) <- "double"

  eta
}

row_max_base <- function(x) {
  x <- as.matrix(x)

  if (ncol(x) == 0L) {
    stop("x must have at least one column.")
  }

  out <- x[, 1]

  if (ncol(x) > 1L) {
    for (j in 2:ncol(x)) {
      out <- pmax(out, x[, j])
    }
  }

  out
}

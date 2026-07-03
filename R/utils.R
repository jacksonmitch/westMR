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

  eta <- (A %*% t(beta_g)) + common_eta(B, beta)
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

format_none <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return("none")
  }

  paste(x, collapse = ", ")
}

coerce_binomial_formula <- function(formula, data) {
  # Just in case
  if (".binom_size" %in% names(data)) {
    stop("Column name '.binom_size' is reserved; please rename your column.")
  }

  resp <- formula[[2]]

  if (is.call(resp) && identical(resp[[1]], as.name("cbind"))) {
    success_expr <- resp[[2]]
    failure_expr <- resp[[3]]

    successes <- eval(success_expr, data)
    failures <- eval(failure_expr, data)

    if (any(successes < 0) || any(failures < 0)) {
      stop("Successes and failures must be non-negative.")
    }

    data$.binom_size <- successes + failures
    formula[[2]] <- success_expr # response becomes the success count only instead of am ugly cbind

    return(list(formula = formula, data = data))
  }

  # Plain response: treat as 0/1 (Bernoulli, size = 1)
  y <- eval(resp, data)

  if (!all(y %in% c(0, 1))) {
    stop("Binomial response must be cbind(successes, failures) or a 0/1 vector.")
  }

  data$.binom_size <- rep(1, length(y))
  list(formula = formula, data = data)
}

compact <- function(x) Filter(Negate(is.null), x)

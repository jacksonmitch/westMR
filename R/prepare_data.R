# Validated data container for mixture regression fitting.
# All structural checks happen once here so callers can trust the contents.

prepare_data <- function(model, included = model$predictors, common = NULL,
                         binomial_size = model$binomial_size) {


  mf <- stats::model.frame(
    formula = model$formula,
    data = model$data,
    na.action = stats::na.fail
  )

  y <- as.numeric(stats::model.response(mf))

  binomial_size_vec <- NULL

  if (identical(model$family, "binomial")) {
    n <- length(y)

    if (is.null(binomial_size)) {
      binomial_size_vec <- rep(1, n)
    } else if (length(binomial_size) == 1L && is.character(binomial_size)) {
      if (!binomial_size %in% names(model$data)) {
        stop("binomial_size variable not found in data.")
      }

      binomial_size_vec <- model$data[[binomial_size]]
    } else if (length(binomial_size) == n) {
      binomial_size_vec <- binomial_size
    }
  binomial_size_vec <- as.numeric(binomial_size_vec)
  }



  het_formula <- make_formula(
    predictors = setdiff(included, common),
    response = model$response
  )
  # Creates intercept column
  X_het <- stats::model.matrix(het_formula, data = mf)

  if (length(common) > 0) {
    com_formula <- stats::as.formula(
      paste("~", paste(common, collapse = " + "), "- 1")
    )
    X_com <- stats::model.matrix(com_formula, data = mf)
  } else {
    X_com <- matrix(numeric(0), nrow = length(y), ncol = 0)
  }

  WMRData$new(y = y, X_het = X_het, X_com = X_com, binomial_size = binomial_size_vec)
}


WMRData <- R6::R6Class(
  "WMRData",
  private = list(
    .y = NULL,
    .X_het = NULL,
    .X_com = NULL,
    .binomial_size = NULL
  ),
  active = list(
    y = function(v) {
      if (missing(v)) private$.y else stop("y is read-only")
    },
    X_het = function(v) {
      if (missing(v)) private$.X_het else stop("X_het is read-only")
    },
    X_com = function(v) {
      if (missing(v)) private$.X_com else stop("X_com is read-only")
    },
    binomial_size = function(v) {
      if (missing(v)) private$.binomial_size else stop("binomial_size is read-only")
    },
    n = function() nrow(private$.X_het),
    p_het = function() ncol(private$.X_het),
    p_com = function() ncol(private$.X_com)
  ),
  public = list(
    initialize = function(y, X_het, X_com, binomial_size = NULL) {
      if (!is.numeric(y)) stop("y must be a numeric vector.")
      if (!is.matrix(X_het)) stop("X_het must be a matrix.")
      if (!is.matrix(X_com)) stop("X_com must be a matrix.")

      n <- length(y)

      if (nrow(X_het) != n) {
        stop(sprintf(
          "nrow(X_het) (%d) must equal length(y) (%d).",
          nrow(X_het), n
        ))
      }
      if (nrow(X_com) != n) {
        stop(sprintf(
          "nrow(X_com) (%d) must equal length(y) (%d).",
          nrow(X_com), n
        ))
      }

      if (n <= ncol(X_het)) {
        stop(sprintf(
          "Need n > ncol(X_het) for the weighted QR M-step.
          Got n = %d, ncol(X_het) = %d.",
          n, ncol(X_het)
        ))
      }

      # na.fail in model.frame covers y and X_het; X_com is built separately
      if (anyNA(X_com)) {
        stop("X_com contains missing values.")
      }
      if (!is.null(binomial_size)) {
        binomial_size <- as.numeric(binomial_size)
      }

      private$.y <- y
      private$.X_het <- X_het
      private$.X_com <- X_com
      private$.binomial_size <- binomial_size
    },
    print = function(...) {
      cat(
        "<WMRData>  n =", self$n,
        " p_het =", self$p_het,
        " p_com =", self$p_com, "\n"
      )
      invisible(self)
    }
  )
)

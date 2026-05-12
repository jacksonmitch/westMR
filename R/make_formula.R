make_rhs <- function(vars, intercept = TRUE) {
  vars <- unique(vars)
  
  if (length(vars) == 0) {
    if (intercept) {
      return("1")
    } else {
      return("0")
    }
  }
  
  rhs <- paste(vars, collapse = " + ")
  
  if (!intercept) {
    rhs <- paste(rhs, "- 1")
  }
  
  rhs
}

make_gmr_formula <- function(response, heterogeneous) {
  stats::as.formula(
    paste(response, "~", make_rhs(heterogeneous, intercept = TRUE))
  )
}

make_common_formula <- function(homogeneous) {
  stats::as.formula(
    paste("~", make_rhs(homogeneous, intercept = FALSE))
  )
}
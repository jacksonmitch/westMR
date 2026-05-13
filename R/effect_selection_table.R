effect_selection_table <- function(object) {
  if (!inherits(object, "effect_selection")) {
    stop("object must be an effect_selection object.")
  }
  
  if (length(object$steps) == 0) {
    return(data.frame())
  }
  
  rows <- lapply(object$steps, function(s) {
    data.frame(
      step = s$step,
      direction = s$direction,
      covariate = names(s$p0),
      p0 = as.numeric(s$p0),
      chosen = names(s$p0) == s$chosen,
      heterogeneous_before = paste(s$heterogeneous_before, collapse = ", "),
      homogeneous_before = paste(s$homogeneous_before, collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, rows)
}
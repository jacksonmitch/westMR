
westMR <- function(
    formula,
    data,
    G_max = 4,
    family = c("gaussian", "poisson"),
    task = c("variables", "effects"),
    control = build_control()
) {

  family <- match.arg(family) # defaults to gaussian

  collection <- checkmate::makeAssertCollection()
  # User Input Checks
  checkmate::assert_formula(formula, add = collection)
  checkmate::assert_data_frame(data, add = collection)
  if (!isTRUE(checkmate::check_int(G_max, lower = 2))) {
    collection$push(sprintf("Maximum number of Components/Groups (G_max) must be >= 2 (received: %s).", G_max))
  }
  checkmate::assert_choice(family, choices = c("gaussian", "poisson"), add = collection)
  checkmate::assert_subset(task, choices = c("variables", "effects"), add = collection)

  checkmate::assert_class(control, "WMRControl", add = collection)
  if (is.null(control$sigma_floor)) {
    mf <- stats::model.frame(
      formula = formula,
      data = data,
    )
    response <- as.numeric(stats::model.response(mf))
    control$sigma_floor <- 0.05 * stats::sd(response)
  }

  if (!collection$isEmpty()) {
    err_messages <- paste0("- ", collection$getMessages(), collapse = "\n")
    stop("(westMR) Error with arguments: \n", err_messages, call. = FALSE)
  }

  # Specifying functions according to distribution family
  if (family == "gaussian"){
    log_lik = obs_loglik_gmr
    m_step_var = m_step_qr
    m_step_eff = m_step_sqr
  }
  else if (family == "poisson"){
    stop("poisson not implemented yet")
  }


  # Create internal model
  model <- WMRModel$new(
    formula = formula,
    data = data,
    G_values = 1:G_max,  # seq_len(as.integer(G_max)) ??
    m_step_var = m_step_var,
    m_step_eff = m_step_eff,
    log_lik = log_lik,
    control = control
  )

  # Run Variable Selection
  if ("variables" %in% task){
    cat("testing variables")
  }
  # Run Effect Type Determination
  if ("effects" %in% task){
    effect_selection <- determine_effects(model, direction = control$direction)
  }
}

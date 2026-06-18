# EM algorithm for finite mixture regression

em_fmr <- function(prepared_data,
                   G,
                   em_state,
                   family,
                   control) {
  # Initial values carried through EM iterations

  loglik_trace <- numeric(control$max_iter)
  converged <- FALSE

  stopifnot(!is.null(em_state[["tau"]]))

  if (is.null(em_state[["pi_g"]])) {
    em_state[["pi_g"]] <- colMeans(em_state[["tau"]])
  }
  irwls_iterations <- integer(control$max_iter)
  irwls_converged <- logical(control$max_iter)

  # Initial log-likelihood before any EM updates
  loglik_old <- -Inf

  for (iter in seq_len(control$max_iter)) {
    # M-step: update parameters using the selected family

    em_state <- m_step(prepared_data, em_state, family, control)

    if (!is.null(em_state[["irwls_iterations"]])) {
      irwls_iterations[iter] <- em_state[["irwls_iterations"]]
    }
    if (!is.null(em_state[["irwls_converged"]])) {
      irwls_converged[iter] <- isTRUE(em_state[["irwls_converged"]])
    }

    # E-step: compute posterior probabilities
    em_state <- e_step_fmr(prepared_data, em_state, family)

    loglik_new <- em_state[["loglik"]]
    loglik_trace[iter] <- loglik_new
    # Convergence check
    tol <- control$tol
    if (is.finite(loglik_old) && tol > 0) {
      ll_diff <- abs(loglik_new - loglik_old)
      ll_scale <- 1 + abs(loglik_old)

      if (ll_diff < tol * ll_scale) {
        converged <- TRUE
        break
      }
    }

    loglik_old <- loglik_new
  }

  loglik_trace <- loglik_trace[seq_len(iter)]
  irwls_iterations <- irwls_iterations[seq_len(iter)]
  irwls_converged <- irwls_converged[seq_len(iter)]

  out <- list(
    em_state = em_state,
    loglik = loglik_trace[length(loglik_trace)],
    loglik_trace = loglik_trace,
    iterations = iter,
    converged = converged,
    irwls_iterations = if (family == "gaussian") NULL else irwls_iterations,
    irwls_converged = if (family == "gaussian") NULL else irwls_converged
  )

  out
}

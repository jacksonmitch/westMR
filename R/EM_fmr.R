# EM algorithm for finite mixture regression

em_fmr <- function(prepared_data, G, em_state, family, control) {
  loglik_trace <- numeric(control$max_iter)
  irwls_iterations <- integer(control$max_iter)
  irwls_converged <- logical(control$max_iter)
  converged <- FALSE
  loglik_old <- -Inf

  for (iter in seq_len(control$max_iter)) {
    m_step(prepared_data, em_state, family, control)

    irwls_iterations[iter] <- em_state$irwls_iterations
    irwls_converged[iter] <- em_state$irwls_converged

    e_step_fmr(prepared_data, em_state, family)

    loglik_new <- em_state$loglik
    loglik_trace[iter] <- loglik_new

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

  list(
    em_state = em_state,
    loglik = loglik_trace[length(loglik_trace)],
    loglik_trace = loglik_trace,
    iterations = iter,
    converged = converged,
    irwls_iterations = if (family == "gaussian") NULL else irwls_iterations,
    irwls_converged = if (family == "gaussian") NULL else irwls_converged
  )
}

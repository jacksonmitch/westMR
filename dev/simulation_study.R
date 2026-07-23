devtools::load_all()

match_groups <- function(tau, true_group) {
  n <- nrow(tau)
  g_hat <- ncol(tau)
  true_g <- max(true_group)

  indicator <- matrix(0, nrow = n, ncol = true_g)
  indicator[cbind(seq_len(n), true_group)] <- 1

  cost <- matrix(0, nrow = g_hat, ncol = true_g)
  for (k in seq_len(g_hat)) {
    for (t in seq_len(true_g)) {
      cost[k, t] <- sum((tau[, k] - indicator[, t])^2)
    }
  }

  if (g_hat <= true_g) {
    match_true <- as.integer(clue::solve_LSAP(cost))          # length g_hat
    match_fitted <- rep(NA_integer_, true_g)
    match_fitted[match_true] <- seq_len(g_hat)
  } else {
    match_fitted <- as.integer(clue::solve_LSAP(t(cost)))     # length true_g
  }

  match_fitted
}

run_replications <- function(params, control = build_control(), verbose = TRUE) {
  results <- vector("list", n_reps)
  
  sim_data <- do.call(simulate_fmr, c(params, list(n = n, seed = seed)))

  for (i in seq_len(n_reps)) {
    if (verbose) cat(sprintf("Replication %d / %d...\n", i, n_reps))
    results[[i]] <- run_one_replication(
      rep_id = i, sim_data = sim_data, params = params, G_max = G_max, 
      control = control
    )
  }
 
  do.call(rbind, results)
}

run_one_replication <- function(rep_id, sim_data, params, G_max, control = build_control()) {
  t0 <- Sys.time()
  fit <- tryCatch(
    westMR(
      formula = sim_data$formula,
      data = sim_data$data,
      G_max = G_max,
      family = "gaussian",
      task = "both",
      control = control
    ),
    error = function(e) e
  )
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
 
  if (inherits(fit, "error")) {
    return(data.frame(
      rep = rep_id, runtime = runtime, error = TRUE,
      error_message = conditionMessage(fit),
      correct_selection = NA, correct_het = NA, correct_hom = NA,
      correct_outcome = NA, G_hat = NA, misclassification_rate = NA,
      rmse_het = NA, rmse_hom = NA
    ))
  }
 
  selected <- fit$variable_selection$selected
  het_hat <- fit$effect_determination$heterogeneous
  hom_hat <- fit$effect_determination$homogeneous
 
  correct_selection <- isTRUE(setequal(selected, union(true_het, true_hom)))
  correct_het <- isTRUE(setequal(het_hat, true_het))
  correct_hom <- isTRUE(setequal(hom_hat, true_hom))
  correct_outcome <- correct_selection && correct_het && correct_hom
 
  true_G <- length(params$pi)
 
  ari <- NA_real_
  rmse_het <- NA_real_
  rmse_hom <- NA_real_
  
  beta_g_hat <- fit$best_fit$parameter_values$beta_g
  beta_hat <- fit$best_fit$parameter_values$beta
  het_truth <- params$betas[, true_het, drop = FALSE]

  G_hat <- fit$best_fit$G
  tau <- fit$best_fit$parameter_values$tau
  predicted_group <- apply(tau, 1, which.max)

  # Tau metrics
  entropy <- -rowSums(tau * log(pmax(tau, 1e-12)))
  normalized_entropy <- mean(entropy) / log(G_hat)
  entropy_R2 <- 1 - normalized_entropy   # closer to 1 = cleaner, more decisive separation
  mean_max_tau <- mean(apply(tau, 1, max))
  
  metrics <- aricode::compare_clustering(predicted_group, sim_data$true_group, AMI = TRUE)
  ari <- metrics$ARI
  ami <- metrics$AMI

  match <- match_groups(tau, sim_data$true_group)
  naive_accuracy <- mean(predicted_group == match[sim_data$true_group])

  indicator <- matrix(0, nrow = nrow(tau), ncol = true_G)
  indicator[cbind(seq_len(nrow(tau)), sim_data$true_group)] <- 1

  tau_aligned <- matrix(0, nrow = nrow(tau), ncol = true_G)
  for (t in seq_len(true_G)) {
    if (!is.na(match[t])) {
      tau_aligned[, t] <- tau[, match[t]]
    }
  }
  brier <- mean(rowSums((tau_aligned - indicator)^2))

  # RMSE's
  errors <- list()

  for (v in true_het) {
    true_vals <- params$betas[, v]
    err <- if (!(v %in% selected)) {
      true_vals - 0
    } else if (v %in% het_hat) {
      fitted_vals <- rep(NA_real_, true_G)
      matched <- !is.na(match)
      fitted_vals[matched] <- beta_g_hat[match[matched], v]
      fitted_vals - true_vals
    } else {
      rep(beta_hat[[v]], true_G) - true_vals
    }
    errors[[v]] <- data.frame(variable = v, true_role = "het", sq_error = err^2)
  }

  for (v in true_hom) {
    true_val <- params$betas[1, v]
    err <- if (!(v %in% selected)) {
      true_val - 0
    } else if (v %in% hom_hat) {
      beta_hat[[v]] - true_val
    } else {
      beta_g_hat[, v] - true_val   # hom mistakenly split into het
    }
    errors[[v]] <- data.frame(variable = v, true_role = "hom", sq_error = err^2)
  }

  errors <- do.call(rbind, errors)

  rmse_het <- sqrt(mean(errors$sq_error[errors$true_role == "het"], na.rm = TRUE))
  rmse_hom <- sqrt(mean(errors$sq_error[errors$true_role == "hom"], na.rm = TRUE))

  data.frame(
    rep = rep_id, runtime = runtime, error = FALSE, error_message = NA,
    correct_selection = correct_selection, correct_het = correct_het,
    correct_hom = correct_hom, correct_outcome = correct_outcome,
    G_hat = G_hat, entropy_R2 = entropy_R2, naive_accuracy = naive_accuracy,
    ari = ari, ami = ami, brier = brier,
    rmse_het = rmse_het, rmse_hom = rmse_hom
  )
}

summarize_results <- function(results) {
  ok <- results[!results$error, ]

  cat(sprintf("Replications: %d  |  Errors: %d\n\n", nrow(results), sum(results$error)))

  cat("--- Runtime ---\n")
  cat(sprintf("Min: %.2fs  |  Mean: %.2fs  |  Median: %.2fs  |  Max: %.2fs\n\n",
    min(ok$runtime), mean(ok$runtime), stats::median(ok$runtime), max(ok$runtime)
  ))

  cat("--- Model Selection ---\n")
  cat(sprintf("Completely Correct outcome:      %.1f%%\n", 100 * mean(ok$correct_outcome, na.rm = TRUE)))
  cat(sprintf("  Correct variable selection: %.1f%%\n", 100 * mean(ok$correct_selection, na.rm = TRUE)))
  cat(sprintf("  Correct heterogeneous set:  %.1f%%\n", 100 * mean(ok$correct_het, na.rm = TRUE)))
  cat(sprintf("  Correct homogeneous set:    %.1f%%\n", 100 * mean(ok$correct_hom, na.rm = TRUE)))
  cat(sprintf("G recovered correctly (G_hat == true G): %.1f%%\n\n",
    100 * mean(ok$G_hat == true_g, na.rm = TRUE)
  ))

  cat("--- Clustering quality (soft & hard, ground-truth-based) ---\n")
  cat("ARI      : chance-corrected agreement between fitted & true groups (1 = perfect, 0 = chance level)\n")
  cat("AMI      : like ARI, but normalized for cluster-size/entropy differences instead of chance alone\n")
  cat("Accuracy : naive proportion correctly labeled (no chance correction)\n")
  cat("Brier    : soft-assignment error, mean squared distance between tau and the true indicator\n")
  cat("         (0 = perfect confident agreement; penalizes uncertain/wrong tau more than hard accuracy would)\n")
  cat(sprintf("  Mean ARI:      %.4f\n", mean(ok$ari, na.rm = TRUE)))
  cat(sprintf("  Mean AMI:      %.4f\n", mean(ok$ami, na.rm = TRUE)))
  cat(sprintf("  Mean Accuracy: %.4f\n", mean(ok$naive_accuracy, na.rm = TRUE)))
  cat(sprintf("  Mean Brier:    %.4f\n\n", mean(ok$brier, na.rm = TRUE)))

  cat("--- Assignment confidence (truth-free) ---\n")
  cat("Entropy R^2: how confidently the algorithm assigns observations to a group\n")
  cat("             (1 = every observation is confidently assigned to one group; 0 = uniform/uninformative)\n")
  cat(sprintf("  Mean entropy R^2: %.4f\n\n", mean(ok$entropy_R2, na.rm = TRUE)))

  cat("--- Coefficient recovery (RMSE) ---\n")
  cat("Misclassifications and estimates compared against the true values\n")
  cat(sprintf("  RMSE heterogeneous coefficients: %.4f\n", mean(ok$rmse_het, na.rm = TRUE)))
  cat(sprintf("  RMSE homogeneous coefficients:   %.4f\n", mean(ok$rmse_hom, na.rm = TRUE)))

  cat("\n--- Split by G recovery ---\n")
  for (metric in c("ari", "ami", "naive_accuracy", "brier", "entropy_R2")) {
    by_g <- tapply(ok[[metric]], ok$G_hat == true_g, mean, na.rm = TRUE)
    cat(sprintf("  %-15s  G correct: %.4f  |  G wrong: %.4f\n",
      metric, by_g["TRUE"], by_g["FALSE"]
    ))
  }
  
  invisible(ok)
}

true_het <- c("het1", "het2", "het3", "het4", "het5")
true_hom <- c("hom1", "hom2", "hom3")
true_null <- c("null1", "null2", "null3", "null4")

g_max <- 7
seed <- 123
n_reps <- 20
n <- 500
params <- scenarios$four_group_twelve_variables
true_g <- length(params$pi)

sim_results <- rbind(sim_results,
   run_replications(params, build_control(parallel = TRUE)))
sim_results
summarize_results(sim_results)


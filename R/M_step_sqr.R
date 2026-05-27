# sQR M-step using the structured QR decomposition

m_step_sqr <- function(A, B, y, tau, sigma_floor = NULL, return_qr_parts = FALSE) {
  A <- as.matrix(A)


  y <- as.numeric(y)
  tau <- as.matrix(tau)

  n <- nrow(A)
  G <- ncol(tau)
  p1 <- ncol(A)
  q <- ncol(B)
  nG <- n * G

  # Update mixing proportions
  pi_g <- colMeans(tau)

  # Stack y and weights in group-by-group order
  y_stack <- rep(y, times = G)
  w_stack <- as.numeric(tau)

  # Storage
  R1_list <- vector("list", G)
  E_list <- vector("list", G)
  b_list <- vector("list", G)

  Bhat <- matrix(0, nG, q)
  Wy <- numeric(nG)

  idx_g <- function(g) {
    ((g - 1) * n + 1):(g * n)
  }

  # Group-wise QR on W_g^(1/2) A
  # Group-wise QR on W_g^(1/2) A
  for (g in seq_len(G)) {
    ii <- idx_g(g)

    yg <- y_stack[ii]
    wg <- w_stack[ii]

    sw <- sqrt(wg)

    Aw <- A * sw
    yw <- yg * sw

    qrA <- qr(Aw)

    R1g <- qr.R(qrA)[seq_len(p1), seq_len(p1), drop = FALSE]

    tmp_y <- qr.qty(qrA, yw)
    b_g <- tmp_y[seq_len(p1)]

    if (q > 0) {
      Bw <- B * sw

      tmp_B <- qr.qty(qrA, Bw)
      E_g <- tmp_B[seq_len(p1), , drop = FALSE]

      Bhat_g <- qr.resid(qrA, Bw)

      Bhat[ii, ] <- Bhat_g
    } else {
      E_g <- matrix(numeric(0), nrow = p1, ncol = 0)
    }

    yhat_g <- qr.resid(qrA, yw)

    R1_list[[g]] <- R1g
    E_list[[g]] <- E_g
    b_list[[g]] <- b_g

    Wy[ii] <- yhat_g
  }

  # QR solve for the common-effect coefficients beta
  if (q > 0) {
    qrB <- qr(Bhat)

    R2 <- qr.R(qrB)[seq_len(q), seq_len(q), drop = FALSE]

    beta <- as.numeric(
      qr.coef(qrB, Wy)[seq_len(q)]
    )
  } else {
    R2 <- matrix(numeric(0), nrow = 0, ncol = 0)
    beta <- numeric(0)
  }

  # Backsolve for group-specific beta_g
  beta_g <- matrix(
    NA_real_,
    nrow = G,
    ncol = p1
  )

  for (g in seq_len(G)) {
    if (q > 0) {
      rhs <- b_list[[g]] - as.vector(E_list[[g]] %*% beta)
    } else {
      rhs <- b_list[[g]]
    }

    beta_g[g, ] <- backsolve(
      R1_list[[g]],
      rhs
    )
  }

  # Update sigma_g after coefficient update
  mu <- sweep(
    A %*% t(beta_g),
    1,
    common_eta(B, beta),
    "+"
  )

  res2 <- (y - mu)^2
  den <- pmax(1e-8, colSums(tau))

  if (is.null(sigma_floor)) {
    sigma_floor <- 0.05 * stats::sd(y)
  }

  sigma_g <- pmax(
    sigma_floor,
    sqrt(colSums(tau * res2) / den)
  )

  out <- list(
    beta = beta,
    beta_g = beta_g,
    pi = pi_g,
    sigma = sigma_g
  )

  if (isTRUE(return_qr_parts)) {
    out$R1_list <- R1_list
    out$E_list <- E_list
    out$R2 <- R2
  }

  out
}

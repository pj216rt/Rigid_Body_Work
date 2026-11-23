# Generate marker data from the rigid body functional mixed model
#
#theta_ik(t) = mu_k(t) + x1_i * beta(t) + b_i(t) + c_ik(t), and markers
# Y = R(w) s + T + off + noise.

# Lab axes (mm): x = lateral (+ toward subject's left), y = anterior,
# z = superior. Body-frame origin = midpoint of the two condylar apices.
# w1 = rotation about the lateral axis (opening/closing; opening is negative),
# w2 = about the anterior axis (roll), w3 = about the vertical axis (lateral swing).

source("rotations.R")

# Pose dimensions (rad, mm), marker points, and tasks
DIM_NAMES <- c("w1", "w2", "w3", "T1", "T2", "T3")
POINTS <- c("incisor", "condyle_L", "condyle_R")
TASKS <- c("OC", "lat_L", "lat_R", "Ant")

# True task mean curves mu_k(t): rows = time, cols = 6 dims
# All tasks share the out-and-back profile sin(pi t); lat_R mirrors lat_L.
true_mu <- function(t, task) {
  g <- sin(pi * t)
  m <- matrix(0, length(t), 6, dimnames = list(NULL, DIM_NAMES))
  switch(task,
    OC     = { m[, "w1"] <- -0.45 * g; m[, "T2"] <- 7 * g;  m[, "T3"] <- -2.5 * g },
    lat_L    = { m[, "w3"] <-  0.12 * g; m[, "T1"] <- 2 * g;  m[, "w1"] <- -0.05 * g; m[, "w2"] <- 0.02 * g },
    lat_R    = { m[, "w3"] <- -0.12 * g; m[, "T1"] <- -2 * g; m[, "w1"] <- -0.05 * g; m[, "w2"] <- -0.02 * g },
    Ant = { m[, "w1"] <- -0.04 * g; m[, "T2"] <- 6 * g;  m[, "T3"] <- -1.5 * g }
  )
  m
}

# True covariate effect curve beta(t) for one standardized covariate x1
# (same for all tasks)
true_beta <- function(t) {
  m <- matrix(0, length(t), 6, dimnames = list(NULL, DIM_NAMES))
  m[, "w1"] <- -0.06 * sin(pi * t)        # higher x1 -> wider opening
  m[, "T2"] <-  1.2 * sin(pi * t)         # more anterior condylar slide
  m[, "w3"] <-  0.02 * sin(2 * pi * t)    # a shape not proportional to the mean
  m
}

# Simulate jaw movement. Returns long marker data (one row per
# subject x task x time x point, plus a rest frame per subject), the
# covariates, the time grid, and the true parameters.
simulate_jaw <- function(n_sub = 20, n_t = 26, sigma = 0.25,
                         sd_b = c(0.04, 0.015, 0.015, 0.8, 0.8, 0.5),
                         sd_c_frac = 0.5, k_re = 5, seed = 1) {
  set.seed(seed)
  t <- seq(0, 1, length.out = n_t)
  Bre <- splines::bs(t, df = k_re, intercept = TRUE)  # basis for random curves

  # 6 x 6 correlation for subject random curves: opening (w1) negatively
  # correlated with anterior slide (T2) and with downward T3
  Cor <- diag(6)
  Cor[1, 5] <- Cor[5, 1] <- -0.6
  Cor[1, 6] <- Cor[6, 1] <-  0.4
  Cor[5, 6] <- Cor[6, 5] <- -0.3
  Sig_b <- diag(sd_b) %*% Cor %*% diag(sd_b)  # covariance of spline coefficients
  Lb <- chol(Sig_b)

  covs <- data.frame(subject = factor(seq_len(n_sub)), x1 = as.numeric(scale(rnorm(n_sub))))
  mu <- lapply(TASKS, function(k) true_mu(t, k)); names(mu) <- TASKS
  beta <- true_beta(t)

  rows <- list()
  for (i in seq_len(n_sub)) {
    # subject geometry in body frame (mm)
    hw <- 50 + rnorm(1, 0, 2)
    s <- rbind(incisor   = c(rnorm(1, 0, 1), 85 + rnorm(1, 0, 3), -35 + rnorm(1, 0, 3)),
               condyle_L = c( hw, rnorm(1, 0, 1), rnorm(1, 0, 1)),
               condyle_R = c(-hw, 0, 0))
    s["condyle_R", 2:3] <- -s["condyle_L", 2:3]   # origin = condylar midpoint
  off <- c(0, -10, 20) + rnorm(3, 0, 2)  # subject's lab offset

  # reference (rest) frame: theta = 0
  Y0 <- s + matrix(off, 3, 3, byrow = TRUE) + matrix(rnorm(9, 0, sigma), 3, 3)
  rows[[length(rows) + 1]] <- data.frame(subject = i, task = "rest", t = 0,
    point = POINTS, x = Y0[, 1], y = Y0[, 2], z = Y0[, 3])

    b_i <- Bre %*% (matrix(rnorm(k_re * 6), k_re, 6) %*% Lb)  # subject curve
    for (k in TASKS) {
      c_ik <- Bre %*% (matrix(rnorm(k_re * 6), k_re, 6) %*% diag(sd_c_frac * sd_b))  # subject-by-task curve
      th <- mu[[k]] + covs$x1[i] * beta + b_i + c_ik
      for (j in seq_len(n_t)) {
        # rotate about the origin, translate, add noise
        R <- expSO3(th[j, 1:3])
        Y <- s %*% t(R) + matrix(th[j, 4:6] + off, 3, 3, byrow = TRUE) +
             matrix(rnorm(9, 0, sigma), 3, 3)
        rows[[length(rows) + 1]] <- data.frame(subject = i, task = k, t = t[j],
                                               point = POINTS, x = Y[, 1], y = Y[, 2], z = Y[, 3])
      }
    }
  }
  dat <- do.call(rbind, rows)
  dat$subject <- factor(dat$subject, levels = levels(covs$subject))
  rownames(dat) <- NULL
  list(dat = dat, covs = covs, t = t,
       truth = list(mu = mu, beta = list(x1 = beta), Sig_b = Sig_b, sigma = sigma))
  }

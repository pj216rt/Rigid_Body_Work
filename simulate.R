 #THE MODEL
#   For subject i, task k, time t, and marker p:
#
#     marker position  =  R(t) %*% s_p  +  T(t)  +  noise
#
#     s_p   = where marker p sits on the jaw, in the jaw's own frame
#             (origin halfway between the two condyles)
#     R(t)  = rotation of the jaw at time t, built from a rotation vector w(t)
#     T(t)  = where the condylar midpoint has moved to at time t
#
#   The 6 motion curves (w1, w2, w3, T1, T2, T3) are what get the covariates:
#
#     motion(t) = task mean curve(t) + x1 * beta(t) + subject random curve(t) + 
# subject-by-task random curve(t)
#
# AXES (mm): x = lateral (positive = subject's left), y = anterior, z = up.
#   w1 = rotation about the left-right axis  (opening is NEGATIVE)
#   w2 = rotation about the front-back axis  (roll)
#   w3 = rotation about the vertical axis    (side-to-side swing)
#   T1, T2, T3 = movement of the condylar midpoint along x, y, z
# =============================================================================
source("rotations.R")
set.seed(2026)

n_subjects <- 20
n_frames   <- 26                       # frames per task
noise_sd   <- 0.25                     # marker noise, mm
times      <- seq(0, 1, length.out = n_frames)
tasks      <- c("OC", "lat_L", "lat_R", "Ant")
markers    <- c("incisor", "condyle_L", "condyle_R")
motion_names <- c("w1", "w2", "w3", "T1", "T2", "T3")

bump <- sin(pi * times) 

#1. True task mean curves (n_frames x 6 matrix per task)
empty <- matrix(0, n_frames, 6, dimnames = list(NULL, motion_names))

task_mean <- list()
task_mean$OC       <- empty
task_mean$OC[, "w1"] <- -0.45 * bump     # about 26 degrees of opening
task_mean$OC[, "T2"] <-  7.0  * bump     # slide 7 mm forward
task_mean$OC[, "T3"] <- -2.5  * bump     

task_mean$lat_L    <- empty
task_mean$lat_L[, "w3"] <- -0.12 * bump
task_mean$lat_L[, "T1"] <-  2.0  * bump
task_mean$lat_L[, "w1"] <- -0.05 * bump

task_mean$lat_R    <- empty
task_mean$lat_R[, "w3"] <-  0.12 * bump
task_mean$lat_R[, "T1"] <- -2.0  * bump
task_mean$lat_R[, "w1"] <- -0.05 * bump

task_mean$Ant      <- empty
task_mean$Ant[, "T2"] <-  6.0  * bump
task_mean$Ant[, "T3"] <- -1.5  * bump
task_mean$Ant[, "w1"] <- -0.04 * bump

#True covariate effect beta(t) for covariate x1 
true_beta <- empty
true_beta[, "w1"] <- -0.06 * bump              # higher x1: wider opening
true_beta[, "T2"] <-  1.2  * bump              # ... and more forward slide
true_beta[, "w3"] <-  0.02 * sin(2 * pi * times)   # a different shape, to be harder

#Random curves
# A random curve = a small B-spline basis times random coefficients.
random_basis <- splines::bs(times, df = 5, intercept = TRUE)   # n_frames x 5

# Subject curves: the 6 dimensions are correlated (6 x 6 covariance).
subject_sd <- c(0.04, 0.015, 0.015, 0.8, 0.8, 0.5)
subject_cor <- diag(6)
subject_cor[1, 5] <- subject_cor[5, 1] <- -0.6    # more opening with more slide
subject_cor[1, 6] <- subject_cor[6, 1] <-  0.4
subject_cor[5, 6] <- subject_cor[6, 5] <- -0.3
subject_cov <- diag(subject_sd) %*% subject_cor %*% diag(subject_sd)
chol_subject_cov <- chol(subject_cov)

# Subject-by-task curves: independent across dimensions, half the size.
curve_sd <- subject_sd / 2

#Generate the markers
covariates <- data.frame(subject = 1:n_subjects,
                         x1 = as.numeric(scale(rnorm(n_subjects))))

rows <- list()
for (i in 1:n_subjects) {

  # Marker positions in the jaw's own frame (rows = markers, cols = x, y, z).
  half_width <- 50 + rnorm(1, 0, 2)
  incisor    <- c(rnorm(1, 0, 1), 85 + rnorm(1, 0, 3), -35 + rnorm(1, 0, 3))
  condyle_L  <- c(half_width, rnorm(1, 0, 1), rnorm(1, 0, 1))
  condyle_R  <- c(-condyle_L[1], -condyle_L[2], -condyle_L[3])   # midpoint = origin
  body <- rbind(incisor, condyle_L, condyle_R)

  # Where this subject's head sits in the lab (an arbitrary offset).
  offset <- c(0, -10, 20) + rnorm(3, 0, 2)

  # Rest frame: no rotation, no translation.
  rest <- body + matrix(offset, 3, 3, byrow = TRUE) + matrix(rnorm(9, 0, noise_sd), 3, 3)
  rows[[length(rows) + 1]] <- data.frame(subject = i, task = "rest", t = 0, marker = markers,
                                         x = rest[, 1], y = rest[, 2], z = rest[, 3])

  # This subject's random curve (n_frames x 6).
  subject_coefs <- matrix(rnorm(5 * 6), 5, 6) %*% chol_subject_cov
  subject_curve <- random_basis %*% subject_coefs

  for (k in tasks) {
    curve_coefs <- matrix(rnorm(5 * 6), 5, 6) %*% diag(curve_sd)
    task_curve  <- random_basis %*% curve_coefs

    # The 6 motion curves for this subject and task.
    motion <- task_mean[[k]] + covariates$x1[i] * true_beta + subject_curve + task_curve

    for (j in 1:n_frames) {
      R <- rotation_vector_to_matrix(motion[j, 1:3])
      T <- motion[j, 4:6]
      lab <- t(R %*% t(body)) +                         # rotate each marker
             matrix(T + offset, 3, 3, byrow = TRUE) +   # then translate
             matrix(rnorm(9, 0, noise_sd), 3, 3)        # then add noise
      rows[[length(rows) + 1]] <- data.frame(subject = i, task = k, t = times[j], marker = markers,
                                             x = lab[, 1], y = lab[, 2], z = lab[, 3])
    }
  }
}
marker_data <- do.call(rbind, rows)

saveRDS(list(markers = marker_data, covariates = covariates,
             truth = list(beta = true_beta, times = times, subject_cor = subject_cor)),
        "sim_data.rds")
cat("Saved sim_data.rds:", nrow(marker_data), "marker rows\n")

# PLOTS of the simulated data
dir.create("plots", showWarnings = FALSE)

#helper function
save_plot <- function(file, draw, width = 1100, height = 900) {
  png(file.path("plots", file), width = width, height = height, res = 110)
  draw()
  dev.off()
  if (interactive()) draw()
  invisible(NULL)
}

task_cols <- c(OC = "black", lat_L = "blue", lat_R = "red", Ant = "darkgreen")

# Displacement of each marker from that subject's rest position.
rest_pos <- marker_data[marker_data$task == "rest", c("subject", "marker", "x", "y", "z")]
names(rest_pos)[3:5] <- c("x0", "y0", "z0")
moving <- merge(marker_data[marker_data$task != "rest", ], rest_pos, by = c("subject", "marker"))
moving$dx <- moving$x - moving$x0
moving$dy <- moving$y - moving$y0
moving$dz <- moving$z - moving$z0
moving <- moving[order(moving$subject, moving$task, moving$marker, moving$t), ]
incisor <- moving[moving$marker == "incisor", ]

# Plot 1: incisor displacement over time, one line per subject.
# Rows = tasks, columns = x (lateral), y (anterior), z (vertical).
save_plot("1_incisor_vs_time.png", function() {
  par(mfrow = c(length(tasks), 3), mar = c(3, 3, 2, 1), mgp = c(2, 0.7, 0))
  for (k in tasks) {
    for (axis in c("dx", "dy", "dz")) {
      d <- incisor[incisor$task == k, ]
      plot(NA, xlim = c(0, 1), ylim = range(incisor[[axis]]),
           xlab = "t", ylab = paste(axis, "(mm)"), main = paste("Incisor", k, axis))
      for (i in 1:n_subjects) {
        di <- d[d$subject == i, ]
        lines(di$t, di[[axis]], col = adjustcolor(task_cols[k], 0.5))
      }
    }
  }
})

# Plot 2: incisor path in space, in the plane where each task moves most.
# OC and Ant: sagittal view (anterior vs vertical). Lateral tasks: top view.
save_plot("2_incisor_paths.png", function() {
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  for (k in tasks) {
    d <- incisor[incisor$task == k, ]
    if (k %in% c("OC", "Ant")) { h <- "dy"; v <- "dz"; hl <- "anterior (mm)"; vl <- "vertical (mm)" }
    else                       { h <- "dx"; v <- "dy"; hl <- "lateral, + = left (mm)"; vl <- "anterior (mm)" }
    plot(NA, xlim = range(d[[h]]), ylim = range(d[[v]]), xlab = hl, ylab = vl,
         main = paste("Incisor path:", k), asp = 1)
    for (i in 1:n_subjects) {
      di <- d[d$subject == i, ]
      lines(di[[h]], di[[v]], col = adjustcolor(task_cols[k], 0.5))
    }
  }
})

# Plot 3: condyle paths (sagittal view), left and right, all tasks.
save_plot("3_condyle_paths.png", width = 1100, height = 600, function() {
  par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  condyles <- moving[moving$marker != "incisor", ]
  for (m in c("condyle_L", "condyle_R")) {
    d <- condyles[condyles$marker == m, ]
    plot(NA, xlim = range(condyles$dy), ylim = range(condyles$dz),
         xlab = "anterior (mm)", ylab = "vertical (mm)", main = paste(m, "path"), asp = 1)
    for (k in tasks) for (i in 1:n_subjects) {
      di <- d[d$task == k & d$subject == i, ]
      lines(di$dy, di$dz, col = adjustcolor(task_cols[k], 0.4))
    }
    legend("bottomleft", legend = tasks, col = task_cols, lwd = 2, bty = "n", cex = 0.8)
  }
})

# Plot 4: the true curves that generated the data:
# task means (solid) and the covariate effect beta (dashed grey).
save_plot("4_true_curves.png", function() {
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  for (d in motion_names) {
    vals <- c(sapply(tasks, function(k) task_mean[[k]][, d]), true_beta[, d])
    plot(NA, xlim = c(0, 1), ylim = range(vals), xlab = "t",
         ylab = if (d %in% c("w1", "w2", "w3")) "rad" else "mm", main = paste("True", d))
    for (k in tasks) lines(times, task_mean[[k]][, d], col = task_cols[k], lwd = 2)
    lines(times, true_beta[, d], col = "grey40", lwd = 2, lty = 2)
    if (d == "w1") legend("bottomleft", legend = c(tasks, "beta"), col = c(task_cols, "grey40"),
                          lty = c(rep(1, 4), 2), lwd = 2, bty = "n", cex = 0.8)
  }
})

# Plot 5: the jaw moving in 3D, one panel per task, for one subject.
# The jaw is drawn as the triangle incisor - condyle_L - condyle_R at a few
# frames (light = start, dark = peak of the movement); thin lines are the
# full paths of the three markers. Positions are relative to the rest frame's
# condylar midpoint, so (0, 0, 0) is where the condyles sit at rest.
show_subject <- 1
save_plot("5_jaw_motion_3d.png", width = 1200, height = 1100, function() {
  s   <- marker_data[marker_data$subject == show_subject, ]
  r0  <- s[s$task == "rest", ]
  ctr <- colMeans(r0[r0$marker != "incisor", c("x", "y", "z")])  # rest condylar midpoint
  s$x <- s$x - ctr[1]; s$y <- s$y - ctr[2]; s$z <- s$z - ctr[3]
  r   <- s[s$task == "rest", ]
  lims <- lapply(s[, c("x", "y", "z")], range)
  snap_frames <- round(seq(1, (n_frames + 1) / 2, length.out = 4))  # start -> peak
  snap_cols   <- colorRampPalette(c("grey80", "black"))(length(snap_frames))
  marker_cols <- c(incisor = "darkorange", condyle_L = "blue", condyle_R = "red")

  par(mfrow = c(2, 2), mar = c(1, 2.5, 2, 1))
  for (k in tasks) {
    # An empty 3D box to draw into (scale = FALSE keeps mm equal on all axes).
    pm <- persp(lims$x, lims$y, matrix(lims$z[1], 2, 2), zlim = lims$z,
                theta = 35, phi = 20, scale = FALSE, border = NA, box = TRUE,
                ticktype = "detailed", cex.axis = 0.6, cex.lab = 0.8,
                xlab = "x lateral (+ left)", ylab = "y anterior", zlab = "z up",
                main = paste0(k, " (subject ", show_subject, ")"))
    to2d <- function(p) trans3d(p$x, p$y, p$z, pm)
    d <- s[s$task == k, ]

    # Full path of each marker over the task.
    for (m in markers) {
      dm <- d[d$marker == m, ]
      lines(to2d(dm[order(dm$t), ]), col = marker_cols[m], lwd = 1.5)
    }
    # The rest-frame jaw (dashed) and snapshots of the moving jaw.
    rr <- r[match(c(markers, markers[1]), r$marker), ]
    lines(to2d(rr), lty = 2, col = "grey50")
    for (f in seq_along(snap_frames)) {
      df <- d[d$t == times[snap_frames[f]], ]
      df <- df[match(c(markers, markers[1]), df$marker), ]
      lines(to2d(df), col = snap_cols[f], lwd = 2)
      points(to2d(df[1:3, ]), pch = 19, cex = 0.8, col = marker_cols[markers])
    }
    if (k == tasks[1])
      legend("bottomright", legend = c(markers, "jaw at rest", "jaw: start -> peak"),
             col = c(marker_cols, "grey50", "black"), lty = c(1, 1, 1, 2, 1), lwd = 2,
             bty = "n", cex = 0.75)
  }
})

cat("Saved plots to", normalizePath("plots"), "\n")

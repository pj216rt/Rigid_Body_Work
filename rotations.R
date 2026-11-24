# The "cross product matrix" of w: a 3 x 3 matrix K with K %*% v = cross(w, v).
cross_matrix <- function(w) {
  rbind(c(    0, -w[3],  w[2]),
        c( w[3],     0, -w[1]),
        c(-w[2],  w[1],     0))
}

# Rotation vector -> rotation matrix (Rodrigues' formula).
rotation_vector_to_matrix <- function(w) {
  angle <- sqrt(sum(w^2))
  K <- cross_matrix(w)
  if (angle < 1e-10) {
    return(diag(3) + K)                  # tiny rotation: R is almost the identity
  }
  diag(3) + (sin(angle) / angle) * K + ((1 - cos(angle)) / angle^2) * (K %*% K)
}

# Rotation matrix -> rotation vector (the inverse of the function above).
rotation_matrix_to_vector <- function(R) {
  cos_angle <- (sum(diag(R)) - 1) / 2
  angle <- acos(cos_angle)

  A <- R - t(R)
  axis_times_2sin <- c(A[3, 2], A[1, 3], A[2, 1])

  if (angle < 1e-8) {
    return(axis_times_2sin / 2)          # tiny rotation
  }
  angle * axis_times_2sin / (2 * sin(angle))
}

# Best rigid fit of one frame (Kabsch algorithm)
fit_rigid_motion <- function(body_points, lab_points) {
  #Center both point sets at their own centroid.
  body_center <- colMeans(body_points)
  lab_center  <- colMeans(lab_points)

  #iterate over matrix and apply summary statistics 
  body_c <- sweep(body_points, 2, body_center)
  lab_c  <- sweep(lab_points,  2, lab_center)

  #SVD of the 3 x 3 cross-covariance.
  H <- t(body_c) %*% lab_c
  s <- svd(H)
  
  #(determinant must be +1).
  flip <- sign(det(s$v %*% t(s$u)))
  R <- s$v %*% diag(c(1, 1, flip)) %*% t(s$u)

  #translation lines up the centroids after rotating.
  T <- lab_center - as.vector(R %*% body_center)

  #report error
  predicted <- t(R %*% t(body_points)) + matrix(T, 3, 3, byrow = TRUE)
  rms_error <- sqrt(mean(rowSums((lab_points - predicted)^2)))

  list(R = R, T = T, rms_error = rms_error)
}
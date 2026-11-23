# SO(3) helpers (rotations as rotation vectors: axis * angle in radians)

# Skew-symmetric cross-product matrix of w
hat3 <- function(w) {
  # [w]x, so that hat3(w) %*% v == cross(w, v)
  matrix(c(0, w[3], -w[2],
           -w[3], 0, w[1],
           w[2], -w[1], 0), 3, 3)
}

# Inverse of hat3
vee3 <- function(A) c(A[3, 2], A[1, 3], A[2, 1])

# Rotation matrix from rotation vector w (Rodrigues' formula);
# first-order approximation near zero
expSO3 <- function(w) {
  th <- sqrt(sum(w^2))
  K <- hat3(w)
  if (th < 1e-10) return(diag(3) + K)
  diag(3) + (sin(th) / th) * K + ((1 - cos(th)) / th^2) * (K %*% K)
}

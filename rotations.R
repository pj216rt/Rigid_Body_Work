# SO(3) helpers

hat3 <- function(w) {
  # [w]x, so that hat3(w) %*% v == cross(w, v)
  matrix(c(0, w[3], -w[2],
           -w[3], 0, w[1],
           w[2], -w[1], 0), 3, 3)
}

vee3 <- function(A) c(A[3, 2], A[1, 3], A[2, 1])
vect <- function(matrix){
  require(lavaan)
  d <- diag(diag(matrix))
  m <- (matrix-diag(diag(matrix)))*sqrt(2)
  lav_matrix_vechr(m+d)
}

RaoDist <- function(matrix1, matrix2){
  (sum(log(eigen(solve(matrix1)%*%matrix2)$values)^2))
}
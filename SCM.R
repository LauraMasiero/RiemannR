SCM<- function(matrix){
  t <- ncol(matrix)
  scm <- (t-1)^(-1)*matrix%*%t(matrix)
  return(scm)
}

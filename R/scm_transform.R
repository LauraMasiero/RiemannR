scm_transform <- function(scm_list, centroid){
  require(powerplus)
  require(expm)
  C_2star <- Matpow(centroid, -0.5,1)
  transformed_matrices <- list()
  cat("Transforming matrices...", "\n")
  for(i in 1:length(scm_list))
    transformed_matrices[[i]] <- logm(C_2star%*%scm_list[[i]]%*%C_2star)
  cat("Vectorizing matrices...", "\n")
  vect_matrices <- lapply(transformed_matrices, vect)
  vect_matrices <- lapply(vect_matrices, as.vector)
  final_matrix<- matrix(NA, length(vect_matrices), length(vect_matrices[[1]]))
  cat("Constructing final matrix...","\n")
  for(i in 1:length(vect_matrices))
    final_matrix[i,] <- unlist(vect_matrices[[i]])[1:length(vect_matrices[[i]])]
  cat("Done.","\n")
  return(final_matrix)
}
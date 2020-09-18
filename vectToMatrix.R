vectToMatrix <- function(scm_vect_list){
  mat <- matrix(NA, length(scm_vect_list), length(scm_vect_list[[1]]))
  for(i in 1:length(scm_vect_list)) mat[i,] <- unlist(scm_vect_list[[i]])
  return(mat)
}
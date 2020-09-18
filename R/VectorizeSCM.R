VectorizeSCM <- function(scm_list){
  require(purrr)
  final_list <- lapply(flatten(scm_list), vect)
  final_list <- lapply(final_list, as.vector)
  return(final_list)
}
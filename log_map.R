log_map <- function(scm_list, C){
  require(powerplus)
  require(expm)
  C_star <- Matpow(C, 0.5,1)
  C_2star <- Matpow(C, -0.5,1)
  mapped_values <- list()
  for(i in 1:length(scm_list))
    mapped_values[[i]] <- C_star%*%logm(C_2star%*%scm_list[[i]]%*%C_2star)%*%C_star
  return(mapped_values)
}
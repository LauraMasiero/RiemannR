exp_map <- function(mapped_values, C){
  require(powerplus)
  require(expm)
  C_star <- Matpow(C, 0.5,1)
  C_2star <- Matpow(C, -0.5,1)
  orig_space_values <- list()
  orig_space_values[[1]] <- C_star%*%expm(C_2star%*%mapped_values%*%C_2star)%*%C_star
  return(orig_space_values)
}
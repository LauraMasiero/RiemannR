extractBlock <- function(scm_st_list){
  dimension <- unlist(lapply(scm_st_list, function(x) dim(x)[1]))
  dim_block <- unlist(lapply(dimension, function(x) x/2))
  blocks <- list()
  for(d in 1:length(dim_block)){
    new_mat <- matrix(NA, dim_block[d], dim_block[d])
    for(i in 1:dim_block[d]){
      for(j in dim_block[d]:dimension[d]){
        new_mat[j-dim_block[d],i] <- scm_st_list[[d]][j,i]
      }
    }
    blocks[[d]] <- new_mat
  }
  return(blocks)
}
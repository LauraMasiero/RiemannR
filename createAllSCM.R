createAllSCM <- function(list, n=8, st=T, tidy=F, col=17) {
  tmp <- list()
  final_list <- list()
  if(tidy==F){
    if(st==T){
      for(j in 1:n){
        tmp <- list()
        for(i in 1:length(list$super_trial[[j]])){
          tmp[[i]] <- SCM(list$super_trial[[j]][[i]])
          final_list[[j]] <- tmp
        }
      }
    }
    else{
      for(j in 1:n){
        tmp <- list()
        for(i in 1:dim(list[[j]]$epochs)[1]){
          tmp[[i]] <- SCM(list[[j]]$epochs[i,-col,])
          final_list[[j]] <- tmp
        }
      }
    }
  }
  else{
    if(st==T){
      for(j in 1:n){
        tmp <- list()
        for(i in 1:length(list$super_trial[[j]])){
          tmp[[i]] <- SCM(list$super_trial[[j]][[i]])
          final_list[[j]] <- tmp
        }
      }
    }
    else{
      for(j in 1:n){
        tmp <- list()
        for(i in 1:dim(list[[j]]$epochs)[1]){
          tmp[[i]] <- SCM(list[[j]]$epochs[i,,])
          final_list[[j]] <- tmp
        }
      }
    }
  }
  return(final_list)
}
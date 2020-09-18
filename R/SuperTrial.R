SuperTrial <- function(list, n=16, length_epochs=129, tidy=F, col=17){
  P <- length(list)
  mean_epoch <- matrix(NA, n, length_epochs)
  super_trial <- list()
  scm_st <- list()
  if(tidy==F){
    for(i in 1:P){
      new_trials <- list()
      indices<- which(list[[i]]$labels==1, arr.ind = T)
      epochs <- list[[i]]$epochs[indices[,2], -col ,]
      N <- sum(list[[i]]$labels)
      for(c in 1:n){
        for(j in 1:129){
          mean_epoch[c,j] <- sum(epochs[,c,j])/N
        }
      }
      for(t in 1:dim(list[[i]]$epochs)[1]){
        new_trials[[t]] <- rbind(mean_epoch,list[[i]]$epochs[t,-col,] )
        super_trial[[i]]<- new_trials
      }
    }
  }
  else{
    for(i in 1:P){
      new_trials <- list()
      indices<- which(list[[i]]$labels==1, arr.ind = T)
      epochs <- list[[i]]$epochs[indices[,2], ,]
      N <- sum(list[[i]]$labels) 
      for(c in 1:n){
        for(j in 1:129){
          mean_epoch[c,j] <- sum(epochs[,c,j])/N
        }
      }
      for(t in 1:dim(list[[i]]$epochs)[1]){
        new_trials[[t]] <- rbind(mean_epoch,list[[i]]$epochs[t,,] )
        super_trial[[i]]<- new_trials
      }
    }
  }
  return(list(mean_epoch=mean_epoch, super_trial=super_trial))
}


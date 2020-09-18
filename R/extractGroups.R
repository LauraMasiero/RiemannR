extractGroups <- function(list, dataset, events=1){
  group_list <- list()
  if(events==1){
    for(t in 1:length(list)) {
      pos <- which(dataset[[t]]$labels==1, arr.ind=T)[,2]
      for(j in 1:length(list[[t]])) group_list[[t]] <- list[[t]][pos]
    }
  }
  else {
    for(t in 1:length(list)) {
      pos <- which(dataset[[t]]$labels==0, arr.ind=T)[,2]
      for(j in 1:length(list[[t]])) group_list[[t]] <- list[[t]][pos]
    }
  }
  return(group_list)
}
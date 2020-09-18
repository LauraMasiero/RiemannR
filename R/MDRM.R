MDRM <- function(scm_list, centroids_list){
  P <- length(scm_list)
  classes <- rep(0,P)
  distances <- array(NA, dim=c(P,length(centroids_list))) #centroids_list= LENGTH 2
  for(m in 1:P){
    for(c in 1:length(centroids_list))
      distances[m,c] <- RaoDist(centroids_list[[c]], scm_list[[m]])
  }
  classes<- apply(distances, 1, which.min)
  class <- rep(0, P)
  for(i in 1:P){
    if(classes[i]==1) class[i] <- 1
  }
  return(list(classes=class, distances=distances))
}

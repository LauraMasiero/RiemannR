LR_classification <- function(mu, sigma, scm_list){
  cat("Constructing lik_events...", "\n")
  lik_events <- lapply(scm_list, function(x) RG_distr(mu[[1]],sigma[1],16, x))
  cat("Constructing lik_0...", "\n")
  lik_0 <- lapply(scm_list, function(x) RG_distr(mu[[2]],sigma[2],16, x))
  lik_ratios <- list()
  classes <- rep(0,length(scm_list))
  for(i in 1:length(scm_list)){
    lik_ratios[[i]] <- log(lik_events[[i]])-log(lik_0[[i]])
    if(lik_ratios[[i]]>=0) classes[i] <- 1
    #cat("Iteration: ", i, "\t")
  }
  return(list(lik_ratios=lik_ratios, classes=classes))
}